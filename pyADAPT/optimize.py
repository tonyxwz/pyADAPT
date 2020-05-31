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
import platform
import multiprocessing as mp
import logging
import logging.handlers
import logging.config
import time

import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import least_squares

from pyADAPT.dataset import DataSet
from pyADAPT.basemodel import BaseModel
from pyADAPT.regularization import default_regularization
from pyADAPT.timeout import TimeOut

ITER = 1
TIMESTEP = 2
OBJFUNC = 3


class Optimizer(object):
    """optimizes an ADAPT model

    naming:
    - self.parameters: the parameters in the model need optimizing
    - self.parameter_trajectory: parameter trajectory in one iteration
    - self.parameter_trajectories_list: list of numpy.ndarray
    - self.parameter_trajectories: xarray.DataArray, trajectories returned
    same for states and fluxes
    """

    def __init__(
        self, model: BaseModel, dataset: DataSet, parameter_names: list, **options
    ):
        self.model = model
        self.dataset = dataset
        self.parameter_names = list(parameter_names)
        self.parameters: pd.DataFrame = self.model.parameters.loc[self.parameter_names]
        self.options = {
            "optimizer": "trf",
            "lambda_r": 1,
            "odesolver": "BDF",
            "sseThres": 1000,  # for natal's init method
            "ss_time": 1000,
            "R": default_regularization,
            "interpolation": "Hermite",
            "verbose": ITER,  # deprecated, use loglevel
            "loglevel": logging.DEBUG,
            "overture_variation": None,  # pascal, natal are options
            "delta_t": 0.1,
            "n_core": mp.cpu_count(),
            "n_iter": 5,
            "seed": 1,
            "initial_parameters": None,  # TODO give alternatives to the model parameter['init']
            "weights": np.ones(len(self.dataset)),
            "timeout": 100,  # seconds
            "max_retry": 100,
        }
        # model don't know which states and fluxes are visible until data is given
        # arrange the states and fluxes in the dataset to be in the same order as the model
        self.dataset.align(self.model.state_order + self.model.flux_order)
        # in van heerden's dataset, glucose(glc) is absent.
        # state_mask will be ['g1p', 'g6p', 'trh', 't6p', 'udg']
        self.state_mask = self.create_mask(self.dataset.names, self.model.state_order)
        self.flux_mask = self.create_mask(self.dataset.names, self.model.flux_order)
        # logger
        self.time_stamp = time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())
        self.options["logging_config_dict"] = {
            "version": 1,
            "formatters": {
                "brief": {
                    "class": "logging.Formatter",
                    "format": "%(processName)-10s: %(message)s",
                },
                "detailed": {
                    "class": "logging.Formatter",
                    "format": "%(asctime)s %(name)-15s %(levelname)-8s %(processName)-18s %(process)-5d %(message)s",
                },
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "formatter": "brief",
                    "level": "INFO",
                },
                "file": {
                    "class": "logging.FileHandler",
                    "filename": f"adapt_{platform.node()}_{self.time_stamp}.log",
                    "mode": "w",
                    "formatter": "detailed",
                    "level": "DEBUG",  # everything goes into the file
                },
                "warnings": {
                    "class": "logging.FileHandler",
                    "filename": f"adapt_warnings_{platform.node()}_{self.time_stamp}.log",
                    "mode": "w",
                    "formatter": "detailed",
                    "level": "WARNING",  # warning and above
                },
            },
            "loggers": {
                "optim": {
                    "handlers": ["file", "console", "warnings"],
                    "level": "DEBUG",  # change this level when necessary
                },
            },
        }
        self.options.update(options)

    def create_mask(self, names_data, names_model):
        mask = list()
        for i in range(len(names_model)):
            mask.append(names_model[i] in names_data)
        return mask

    def set_timestep(self, delta_t):
        # the end time might not be the actual end time
        self.time = np.arange(self.dataset.begin_time, self.dataset.end_time, delta_t)
        self.options["n_ts"] = len(self.time)

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
        logging.config.dictConfig(self.options["logging_config_dict"])
        logger = logging.getLogger("optim")
        self.set_timestep(self.options["delta_t"])

        logger.info("optimize started")
        logger.info("n_ts:%d", self.options["n_ts"])
        logger.info("n_iter:%d", self.options["n_iter"])
        if self.options["initial_parameters"] is not None:
            self.parameters.loc[:, "init"] = self.options["initial_parameters"]

        self.parameter_trajectories_list = []
        self.state_trajectories_list = []
        self.flux_trajectories_list = []

        man = mp.Manager()
        self.q = man.Queue()

        if self.options["n_core"] > 1:
            pool = mp.Pool(self.options["n_core"])
            logger_thread = threading.Thread(target=self.logger_thread)
            pool_results = []
            logger_thread.start()  # start logging thread
            for i_iter in range(self.options["n_iter"]):
                pool_results.append(
                    pool.apply_async(self.fit_iteration, kwds={"i_iter": i_iter})
                )
            pool.close()
            pool.join()
            self.q.put(None)  # stop logger thread
            logger_thread.join()

            for res_obj in pool_results:
                ptraj, straj, vtraj = res_obj.get()
                self.parameter_trajectories_list.append(ptraj)
                self.state_trajectories_list.append(straj)
                self.flux_trajectories_list.append(vtraj)

        else:
            logger.warning("running large models on single process can be very slow")
            for i_iter in range(self.options["n_iter"]):
                ptraj, straj, vtraj = self.fit_iteration(i_iter=i_iter, parallel=False)
                self.parameter_trajectories_list.append(ptraj)
                self.state_trajectories_list.append(straj)
                self.flux_trajectories_list.append(vtraj)

        # TODO define a class for trajectory
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
        logger.info("optimizer quit successfully")
        return (
            self.parameter_trajectories,
            self.state_trajectories,
            self.flux_trajectories,
        )

    def create_empty_trajectories(self, n_ts):
        self.parameter_trajectory = np.zeros((n_ts, len(self.parameter_names)))
        self.state_trajectory = np.zeros((n_ts, len(self.model.state_order)))
        self.flux_trajectory = np.zeros((n_ts, len(self.model.flux_order)))

    def fit_iteration(self, i_iter=0, parallel=True):
        """ handle one iteration (one subprocess)
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
        # maybe a better way to guaranttee unique seed per process
        np.random.seed(self.options["seed"] + i_iter * 200)
        logger.info("iteration: %d", i_iter)

        for attempt in range(self.options["max_retry"]):
            try:
                self.create_empty_trajectories(self.options["n_ts"])
                splines = self.dataset.interpolate(n_ts=self.options["n_ts"])
                self.overture(i_iter, splines)

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
                break  # break if not a single timestep timeouts
            except TimeoutError as toerr:
                logger.warning(
                    "iter %d timeouts at ts %s on attempt %d",
                    i_iter,
                    str(toerr),
                    attempt,
                )
        return (self.parameter_trajectory, self.state_trajectory, self.flux_trajectory)

    def overture(self, i_iter, splines):
        """optimization before optimization, **init fit**, selling point in this effort
        """
        logger = logging.getLogger("optim.iter.overture")
        with TimeOut(self.options["timeout"], "0") as _timeout:
            # splines: [ observable states | observable fluxes ]
            self.state_trajectory[0, self.state_mask] = splines[
                : sum(self.state_mask), 0, 0
            ]
            for i, observable in enumerate(self.state_mask):
                if not observable:
                    # unobservable state trajectories are left for ADAPT completely
                    self.state_trajectory[0, i] = (
                        self.model.initial_states[i] * np.random.rand()
                    )
            param_init = self.parameters.loc[:, "init"]
            # print("param_init:", param_init)
            op = least_squares(
                self.overture_obj,
                param_init,
                args=(splines[: sum(self.state_mask), 0, 0],),
                bounds=(self.parameters["lb"], self.parameters["ub"]),
                method=self.options["optimizer"],
            )
            if op.success:
                level = logging.DEBUG
            else:
                level = logging.WARNING
            self.parameter_trajectory[0, :] = op.x
            self.parameters.loc[:, "value"] = op.x
            if self.model.has_flux:
                self.flux_trajectory[0, :] = self.model.fluxes(
                    0, self.state_trajectory[0, :], self.model.parameters["value"]
                )
            logger.log(level, "init iter %d %s", i_iter, op.message)
        return op

    def overture_obj(self, params, target):
        """objective function at initial time step

        find **parameters** which
        1. minimize the errors between sampled states and the (long term) steady
            states of the model.
        2. pushing dy to 0 (steady states)

        there are two proposed ways to do this:
        1. initial value problem.
        (might not guarantee the steady states are seamlessly connected to the splines)
            1. solve the ODE on interval [-(ss_time+1), -1] -> sol
            2. get sol.y at t = -1
            3. calculate the fluxes at -1
            4. return residuals
        2. boundary value problem (dropped)
        (looks promising but I haven't learned to use the API, and it doesn't guarantee
         solution existance, at least analytically)
            solve_ivp on [-(sstime+1), -1]
                left boundary: f(-sstime) = A
                right boundary: f(-1) = B,
                                f'(-1) = 0 (fluxes)
        Limitation: there's no support for fluxes for now
        """
        # sstime: how long to simulate the steady state process
        # taking the first approach
        y = self.model.compute_states(
            new_params=params,
            new_param_names=self.parameter_names,
            x0=self.model.initial_states,
            time_points=[self.time[0] - self.options["ss_time"], self.time[0]],
            odesolver=self.options["odesolver"],
        )[:, -1]
        dy = self.model.state_ode(self.time[0], y, self.model.parameters["value"])
        # TODO add extra weighting to the concatenation
        weight = 3
        return np.r_[y[self.state_mask] - target, weight * dy]

    def overture_variation_lazy(self, i_iter, splines):
        """This variation only does no optimzation at all, simply put the random value
        in the trajectory. This is useful when applying the "try many times and select
        only the best x results" strategy.
        """
        logger = logging.getLogger("optim.iter.init")
        with TimeOut(self.options["timeout"], "0") as _timeout:
            # splines: |observable states | observable fluxes|
            self.state_trajectory[0, self.state_mask] = splines[
                : sum(self.state_mask), 0, 0
            ]
            for i, observable in enumerate(self.state_mask):
                if not observable:
                    self.state_trajectory[0, i] = (
                        self.model.initial_states[i] * np.random.rand()
                    )
            # This is risky but according 3-σ principle, this is "almost" always safe
            param_init = np.random.normal(
                self.parameters["init"], self.parameters["init"] * 0.2
            )
            # need to update both because the calculation uses the dataframe while
            # returning trajectory needs the ndarray
            self.parameter_trajectory[0, :] = param_init
            self.parameters.loc[:, "value"] = param_init

            if self.model.has_flux:
                self.flux_trajectory[0, :] = self.model.fluxes(
                    0, self.state_trajectory[0, :], self.model.parameters["value"]
                )
            logger.debug("init iter %d", i_iter)

    def overture_var_pascal(self, i_iter, splines):
        """ Overture variation Pascal
        The initial parameters and states finding method as described in
        Pascal van Beek's master thesis in 2018:
            "Expanding ADAPT to model heterogeneous datasets and application
                to hyperinsulinemic euglycemic clamp data"
        The difference between this method and my method in `initialize_trajectories`
        is that here the parameters are optimized on the entire time span while in my
        method, the optimization happens before the time span.
        """
        _logger = logging.getLogger("optim.iter.overture")

    def overture_var_natal(self, i_iter, splines):
        """ Overture variation Natal
        Basically, it tries multiple times until find a set of parameters that
        fits the data well.
        return: arrays of initial parameters
        """
        _logger = logging.getLogger("optim.iter.overture")
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
                method=self.options["optimizer"],
                kwargs={
                    "begin_states": splines[:, i_ts, 0],
                    "interp_data": splines[:, i_ts, :],
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
        with TimeOut(self.options["timeout"], str(i_ts)) as timeout:
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
                method=self.options["optimizer"],
                **lsq_options,
            )
        # calculate the states and fluxes at the end of the time step using the
        # optimized parameters
        final_states = self.model.compute_states(
            new_params=lsq_result.x,
            time_points=time_span,
            x0=begin_states,
            new_param_names=self.parameter_names,
            odesolver=self.options["odesolver"],
        )[:, -1]
        if self.model.has_flux:
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
            ==================
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

        end_states = self.model.compute_states(
            new_params=params,
            time_points=time_span,
            x0=begin_states,
            new_param_names=self.parameter_names,
            odesolver=self.options["odesolver"],
        )
        end_states = end_states[:, -1]
        if len(self.model.flux_order):
            end_fluxes = self.model.fluxes(
                time_span[-1], end_states, self.model.parameters["value"]
            )
            end_fluxes = end_fluxes[self.flux_mask]
        end_pred = end_states[self.state_mask]

        if len(self.model.flux_order):
            end_pred = np.r_[end_pred, end_fluxes]

        # equation 3.4 [ADAPT 2013]
        residual = (end_pred - end_data[:, 0]) / end_data[:, 1]
        residual = residual * self.options["weights"]
        if self.options["R"] is not None:
            reg_term = self.options["lambda_r"] * self.options["R"](
                params=params,
                parameter_trajectory=self.parameter_trajectory,
                state_trajectory=self.state_trajectory,
                i_iter=i_iter,
                i_ts=i_ts,
                time_span=time_span,
            )
            residual = np.r_[residual, reg_term]
        return residual

    def test_helper(self):
        pass


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
    import pyADAPT.visualize as vis
    import matplotlib.pyplot as plt

    plt.style.use(["science", "grid"])

    model = ToyModel()
    data = DataSet(
        raw_data_path="data/toyModel/toyData.mat",
        data_specs_path="data/toyModel/toyData.yaml",
    )
    ptraj, straj, vtraj, time = optimize(
        model,
        data,
        "k1",
        n_iter=10,
        delta_t=0.2,
        odesolver="LSODA",
        n_core=4,
        verbose=ITER,
        weights=np.array([1, 0.1, 1, 0.5]),
    )

    fig, axes = plt.subplots(figsize=(10, 8))
    vis.plot(ptraj, axes=axes, color="green", alpha=0.2)
    vis.plot_mean(ptraj, axes=axes, color="red")

    # # plt.show()
    # fig.savefig("toy.pdf")
    fig.savefig("test_optimize.png")
