# -*- coding: utf-8 -*-
"""
optimize module just runs ADAPT [^1]
===============================

Motivation
----------

Our system of interest is composed of both the information from metabolic
network and the gene and protein networks. Although metabolic networks are well
studied, a lot of information is still missing from gene and protein networks,
including the *kinetic information*. In order to simulate the system without
the kinetic info of genome, transcriptome and proteome, time-dependent
parameters are introduced to the system. The changes in the levels that we
don't understand such as the proteome and transcriptome is summarized into the
change in time-dependent parameters. The task then becomes finding a parameter/
(parameters)'s corresponding function θ(t), that can best fit the experimental
data.

Part of the reason why ADAPT is so confusing is that neither the paper nor the
various attempt of implementation gives a clear problem solving step of the
algorithm. The paper is ambiguous and the MATLAB® implementations are only
tailored for the authors' own need. I hope this project can change this mess.

What do we have?
----------------

1. experimental data from a true system
2. kinetic information of the metabolic network (incomplete and possibly flawed)

We would like to use the time dependent parameter to compensate the possible
errors in our model.

Procedure
---------

Foundamentally, this is a minimization problem, but separated into small time
steps. The first thing we need to consider is still to define the objective
function ρ. Following the guideline of `scipy.optimize.least_squares` and
`lmfit.minimize` (https://lmfit.github.io/lmfit-py/intro.html), ρ should have a
certain call signature and return value.

- ρ (params, *args, **kws)
- return: residual array
- the steady state of the system as a initial point (t==0)
- the initial parameters θ₀, the initial parameters should be encoded by the
    data provider into the dataset, which guarantees that the first record of
    the data set (at t=0) is the steady state simulation of the project. If
    there's no such data set, then we have to fake some. Fortunately, the toy
    and trehalose model meets such requirements. The first time step data is
    the steady state of the model simulated using the initial parameters.

Program workflow
----------------

optimize
    optim = new Optimizer(model, data, parameter_names)
        optim.run(iteration_number)
            for i in iteration_number:
                for ts in range(n_ts):
                    optim.fit_timestep()
                        scipy:least_square()
                            optim.objective_function
                                optim.model.compute_states
                                    scipy:solve_ivp()
print & plotting

[^1]: I am considering renaming ADAPT to Optimizer in order to be distant away
from this lame paper.
"""

import threading
import multiprocessing as mp
import logging
import logging.handlers
import logging.config
import time

import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import least_squares, leastsq

from pyADAPT.dataset import DataSet
from pyADAPT.basemodel import BaseModel
from pyADAPT.regularization import default_regularization


ITER = 1
TIMESTEP = 2
OBJFUNC = 3


class Optimizer(object):
    """optimizes an ADAPT model

    naming:
    ------

    - self.parameters: the parameters in the model need optimizing
    - self.parameter_trajectory: parameter trajectory in one iteration
    - self.parameter_trajectories_list: list of numpy.ndarray
    - self.parameter_trajectories: xarray.DataArray, trajectories returned

    same for states and fluxes
    """

    def __init__(self, model: BaseModel, dataset: DataSet, parameter_names: list):
        self.model = model
        self.dataset = dataset
        self.parameter_names = list(parameter_names)
        self.parameters: pd.DataFrame = self.model.parameters.loc[self.parameter_names]
        self.options = {
            "method": "trf",
            "lambda_r": 1,
            "odesolver": "RK45",
            "sseThres": 1000,
            "ss_time": 1000,
            "R": default_regularization,
            "interpolation": "Hermite",
            "verbose": ITER,
            "loglevel": logging.DEBUG,
            "init_method": None,  # pascal, natal are options
            "delta_t": 0.1,
            "n_core": mp.cpu_count(),
            "n_iter": 5,
            "seed": 1,
            "weights": np.ones(len(self.dataset)),  # TODO weight of the errors
        }

        # model don't know which states and fluxes are visible until data is given
        # 1. and arrange the states and fluxes in the model to be in the same order as the dataset
        self.dataset.align(self.model.state_order + self.model.flux_order)
        # 2. create the mask here
        #! in van heerden's dataset, glucose(glc) is absent.
        # state_mask will be ['g1p', 'g6p', 'trh', 't6p', 'udg']
        self.state_mask = self.create_mask(self.dataset.names, self.model.state_order)
        self.flux_mask = self.create_mask(self.dataset.names, self.model.flux_order)

        # logger
        self.time_stamp = time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())

        self.options["logging_config_dict"] = {
            "version": 1,
            "formatters": {
                'brief': {
                    'class': 'logging.Formatter',
                    'format': '%(processName)-10s %(message)s'
                },
                'detailed': {
                    'class': 'logging.Formatter',
                    'format': '%(asctime)s %(name)-15s %(levelname)-8s %(processName)-10s %(message)s'
                }
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "formatter": 'brief',
                    "level": "INFO"
                },
                "file": {
                    "class": "logging.FileHandler",
                    "filename": f"adapt_{self.time_stamp}.log",
                    "mode": "w",
                    "formatter": "detailed",
                    "level": "DEBUG"
                }
            },
            "loggers": {
                "optim": {
                    "handlers": ["file", "console"],
                    "level": "DEBUG"   # change this level when necessary
                },
            }
        }

    def create_mask(self, names_data, names_model):
        mask = list()
        for i in range(len(names_model)):
            mask.append(names_model[i] in names_data)
        return mask

    def logger_thread(self):
        while True:
            record = self.q.get()
            if record is None:
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)

    def run(self, **options):
        """ Main Process
        """
        self.options.update(options)
        self.time = np.arange(
            self.dataset.begin_time, self.dataset.end_time, self.options["delta_t"]
        )
        self.options["n_ts"] = len(self.time)
        logging.config.dictConfig(self.options["logging_config_dict"])
        logger = logging.getLogger("optim")
        self.parameter_trajectories_list = []
        self.state_trajectories_list = []
        self.flux_trajectories_list = []

        # start optimizing (in parallel)
        if self.options["n_core"] > 1:
            man = mp.Manager()
            self.q = man.Queue()
            pool = mp.Pool(self.options["n_core"])

            logger_thread = threading.Thread(target=self.logger_thread)
            pool_results = []
            for i_iter in range(self.options["n_iter"]):
                pool_results.append(
                    # make use of the main process's async apply waiting time
                    # and receive logs from other subprocesses
                    pool.apply_async(self.fit_iteration, kwds={"i_iter": i_iter})
                )
            logger_thread.start()  # start logging thread
            pool.close()
            pool.join()
            self.q.put(None)  # break out from logging loop
            logger_thread.join()

            for res_obj in pool_results:
                ptraj, straj, vtraj = res_obj.get()
                self.parameter_trajectories_list.append(ptraj)
                self.state_trajectories_list.append(straj)
                self.flux_trajectories_list.append(vtraj)

        else:
            # TODO just drop the support for single processing, it would be much easier
            for i_iter in range(self.options["n_iter"]):
                ptraj, straj, vtraj = self.fit_iteration(i_iter=i_iter, parallel=False)
                self.parameter_trajectories_list.append(ptraj)
                self.state_trajectories_list.append(straj)
                self.flux_trajectories_list.append(vtraj)

        # convert the lists into xarray
        self.parameter_trajectories = xr.DataArray(
            data=np.array(self.parameter_trajectories_list),
            coords=[
                ("iter", list(range(self.options["n_iter"]))),
                ("time", self.time),
                ("id", list(self.parameter_names)),
            ],
            name="parameter trajectories",
        )
        self.state_trajectories = xr.DataArray(
            data=np.array(self.state_trajectories_list),
            coords=[
                ("iter", list(range(self.options["n_iter"]))),
                ("time", self.time),
                ("id", list(self.model.state_order)),
            ],
            name="state trajectories",
        )
        self.flux_trajectories = xr.DataArray(
            data=np.array(self.flux_trajectories_list),
            coords=[
                ("iter", list(range(self.options["n_iter"]))),
                ("time", self.time),
                ("id", list(self.model.flux_order)),
            ],
            name="flux trajectories",
        )

        return (
            self.parameter_trajectories,
            self.state_trajectories,
            self.flux_trajectories,
        )

    def fit_iteration(self, i_iter=0, parallel=True):
        """ handle one iteration (one subprocess)
        """

        """
        TODO I could write a chapter in the final thesis about how logging works.

        `logger`s in python are singletons, in other words, there's always only
        one logger instance with the same name in one process, and logging.getLogger
        always returns the reference to the same logger object as long as the name
        is the same.

        Things are a bit different in multiprocessing, since loggers are used
        in different processes, different logger object with the same name can
        live inside different processes. The logging system in pyADAPT supports
        multiprocessing by launching a separate thread in the main process to
        handle the log records. Because there are many iterations in one simulation
        and each simulation is processed in a process. The parallelism in pyADAPT
        employs process pool to finish the jobs. As mentioned before, loggers are
        singletons. When a worker process is running a second time, there's no
        need to add handlers to the logger anymore. Or the log will contain
        many duplicated records.
        """

        qh = logging.handlers.QueueHandler(self.q)
        logger = logging.getLogger("optim.iter")
        if not logger.hasHandlers():
            logger.setLevel(logging.DEBUG)
            logger.addHandler(qh)

        #! reseed the random generator of subprocess
        np.random.seed(self.options["seed"] + i_iter)
        # logger.info("reseeded rng")

        logger.info(f"iteration: %d", i_iter)

        self.parameter_trajectory = np.zeros(
            (self.options["n_ts"], len(self.parameter_names))
        )
        self.state_trajectory = np.zeros(
            (self.options["n_ts"], len(self.model.state_order))
        )
        self.flux_trajectory = np.zeros(
            (self.options["n_ts"], len(self.model.flux_order))
        )

        splines = self.dataset.interpolate(n_ts=self.options["n_ts"])
        self.init_ts0(i_iter, splines)

        for i_ts in range(1, self.options["n_ts"]):
            (
                self.parameter_trajectory[i_ts, :],
                self.state_trajectory[i_ts, :],
                self.flux_trajectory[i_ts, :],
            ) = self.fit_timestep(
                initial_guess=self.parameter_trajectory[i_ts - 1, :],
                begin_states=self.state_trajectory[i_ts - 1, :],
                end_data=splines[:, i_ts, :],
                i_iter=i_iter,
                i_ts=i_ts,
            )

        return (self.parameter_trajectory, self.state_trajectory, self.flux_trajectory)

    def init_ts0(self, i_iter, splines):
        """ time step 0 is handled differently because there's no previous
        time step information at this moment

        what need to be done at t=0?
        1. generate state
        2. generate parameter
        3. calculate fluxes from new parameters and states
        """
        logger = logging.getLogger("optim.iter.init")
        self.state_trajectory[0, self.state_mask] = splines[: len(self.model.state_order), 0, 0]
        for i, observable in enumerate(self.state_mask):
            if not observable:
                self.state_trajectory[0, i] = self.model.initial_states[i] * np.random.rand()
        logger.debug("iter %d initialized", i_iter)

        # This is risky but according 3-σ principle, this is "almost" always safe
        param_init = np.random.normal(
            self.parameters["init"], self.parameters["init"] * 0.2
        )
        self.parameter_trajectory[0, :] = param_init
        self.parameters.loc[:, "value"] = param_init
        # TODO add support for self.options['init']
        if self.model.flux_order:
            self.flux_trajectory[0, :] = self.model.fluxes(
                0, self.state_trajectory[0, :], self.model.parameters["value"]
            )

    def pascal_init(self, i_iter, data):
        """ The initial parameters and states finding method inspired by
        Pascal van Beek's master thesis in 2018:
            "Expanding ADAPT to model heterogeneous datasets and application
                to hyperinsulinemic euglycemic clamp data"
        """
        logger = logging.getLogger("optim.iter.pascal")

    def natal_init(self, i_iter, data):
        """ Natal van Riel's parameter initialization
        Basically, it tries multiple times to fit the parameters to the data
        return: arrays of initial parameters
        """
        logger = logging.getLogger("optim.iter.natal")
        sse = np.inf
        i_ts = 0

        while sse > self.options["sseThres"]:
            params = self.parameters["init"] * (
                10
                ** (2 * np.random.random_sample(size=(len(self.parameter_names),)) - 1)
            )
            lsq_res = least_squares(
                self.objective_function,
                params,
                bounds=(self.parameters["lb"], self.parameters["ub"]),
                method=self.options["method"],
                kwargs={
                    "begin_states": data[:, i_ts, 0],
                    "interp_data": data[:, i_ts, :],
                    "time_span": [self.time[0] - self.options["ss_time"], self.time[0]],
                    "R": None,
                    "i_iter": i_iter,
                    "i_ts": i_ts,
                },
            )
            sse = lsq_res.cost
        return lsq_res.x

    def fit_timestep(
        self,
        initial_guess=None,
        begin_states=None,
        end_data=None,
        i_iter=None,
        i_ts=None,
        **lsq_options,
    ):
        """call least_squares
        access Optimizer options via `i_iter` and `i_ts` and `self`
        """
        logger = logging.getLogger("optim.iter.ts")
        logger.debug("iter: %d, ts: %d", i_iter, i_ts)

        time_span = self.time[i_ts - 1 : i_ts + 1]  # just two points
        lsq_result = least_squares(
            self.objective_function,
            initial_guess,
            kwargs={
                "begin_states": begin_states,
                "end_data": end_data,
                "time_span": time_span,
                "i_iter": i_iter,
                "i_ts": i_ts,
            },
            bounds=(self.parameters["lb"], self.parameters["ub"]),
            method=self.options["method"],
            **lsq_options,
        )
        # calculate the states and fluxes at the end of the time step using the
        # optimized parameters
        final_states = self.model.compute_states(
            new_params=lsq_result.x,
            time_points=time_span,
            x0=begin_states,
            new_param_names=self.parameter_names,
        )[:, -1]
        if self.model.flux_order:
            final_fluxes = self.model.fluxes(
                time_span[-1], final_states, self.model.parameters["value"]
            )
        else:
            final_fluxes = []
        return (lsq_result.x, final_states, final_fluxes)

    def objective_function(
        self,
        params,
        begin_states=None,
        end_data=None,
        time_span=None,
        i_iter=None,
        i_ts=None,
        **kw,
    ):
        """ Objective function

            The function minimized by least squares method. For ADAPT, an objective
            should:

            1. compute the states at the end of the `timespan`, using the give
                `parameters`, `begin_states`.
            2. choose those `end_states` and `interp_states` which are "observable"
            3. calculate and choose observable fluxes
            4. calculate residual
            5. calculate regularization term
            6. concatenate errors and regularization terms

            Parameters
            ----------
            params
                the parameters θ of the model, or x in scipy's term
                it should be the length of the parameters to be varied
            begin_states
                the begin states of this time step (end of previous time step)
            interp_data
                shape(len(states), 2), states, stds
                the interpolated experimental DATA at the end of the time span,
                containing the states and corresponding standard deviations.
            time_span
                the time span to solve the ODE
            parameter_names
                list of strings, names of the parameters to optimize, other parameters
                are set as fixed (constant) parameters

            Return
            ------
            np.ndarray: shape(len(states) + len(parameter penalty))
        """
        logger = logging.getLogger("optim.iter.objective")
        if self.options["verbose"] >= OBJFUNC:
            pass  # print something about optimizer

        end_states = self.model.compute_states(
            new_params=params,
            time_points=time_span,
            x0=begin_states,
            new_param_names=self.parameter_names,
        )
        end_states = end_states[:, -1]
        if len(self.model.flux_order):
            end_fluxes = self.model.fluxes(
                time_span[-1], end_states, self.model.parameters["value"]
            )
            end_fluxes = end_fluxes[self.flux_mask]
        end_pred = end_states[self.state_mask]

        if len(self.model.flux_order):
            end_pred = np.concatenate([end_pred, end_fluxes])

        # equation 3.4 [ADAPT 2013]
        residual = (end_pred - end_data[:, 0]) / end_data[:, 1]
        if self.options["R"] is not None:
            reg_term = self.options["lambda_r"] * self.options["R"](
                params=params,
                parameter_trajectory=self.parameter_trajectory,
                state_trajectory=self.state_trajectory,
                i_iter=i_iter,
                i_ts=i_ts,
                time_span=time_span,
            )
            residual = np.concatenate([residual, reg_term])
        return residual


def optimize(model, dataset, *params, **options):
    """ the main optimization procedure

    Parameter
    ---------
    model: a model instance of subclass of BaseModel
    dataset: dataset
    params: list of parameters to be optimized
    n_core: number of processes to spawn

    Return
    ------
    (parameter trajectories, state trajectories, time)
    trajectories are xarray.DataArray. just easier to manipulate in the analysis
    routine.
    """
    optim = Optimizer(model, dataset, params)
    optim.run(**options),

    return (
        optim.parameter_trajectories,
        optim.state_trajectories,
        optim.flux_trajectories,
        optim.time,
    )


if __name__ == "__main__":
    from pyADAPT.examples.lotka import LotkaVolterra
    from pyADAPT.examples.toy import ToyModel
    from pyADAPT.dataset import DataSet
    import pyADAPT.trajectory as traj
    import matplotlib.pyplot as plt

    model = ToyModel()
    data = DataSet(
        raw_data_path="data/toyModel/toyData.mat",
        data_specs_path="data/toyModel/toyData.yaml",
    )
    ptraj, straj, vtraj, time = optimize(
        model, data, "k1", n_iter=10, delta_t=0.2, n_core=4, verbose=ITER
    )

    fig, axes = plt.subplots()
    traj.plot(ptraj, axes=axes, color="green", alpha=0.2)
    traj.plot_mean(ptraj, axes=axes, color="red")

    # plt.show()
