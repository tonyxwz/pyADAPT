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

import datetime
import multiprocessing as mp
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import least_squares, leastsq

from pyADAPT.dataset import DataSet
from pyADAPT.basemodel import BaseModel


def default_regularization(params=None,
                           parameter_trajectory=None,
                           time_span=None,
                           i_iter=None,
                           i_ts=None,
                           **kw):
    """ tiemann & natal's regularization term in ADAPT 2013 paper """
    old_params = parameter_trajectory[i_ts - 1, :]
    delta_t = time_span[-1] - time_span[0]
    reg = (params - old_params) / delta_t / old_params
    return reg


class ADAPTResult(object):
    def __init__(self):
        pass


class Optimizer(object):
    """optimizes an ADAPT model"""
    ITER = 1
    TIMESTEP = 2
    OBJFUNC = 3

    def __init__(self, model: BaseModel, dataset: DataSet,
                 parameter_names: list):
        # I am being naive by assuming the user will give the states in the
        # model and the dataset in the same order, maybe warn in the manual.
        self.model = model
        self.dataset = dataset
        self.parameter_names = parameter_names
        self.parameters: pd.DataFrame = self.model.parameters.loc[list(
            parameter_names)]
        self.options = {  # TODO: a method to set options
            "method": "trf",
            "lamda": 1,
            "odesolver": "RK45",
            "sseThres": 1000,
            "regularization": default_regularization,
            "interpolation": "Hermite",
            "verbose": self.ITER,
            "init_method": None
        }

    def run_mp(self,
               begin_time=None,
               end_time=None,
               n_iter=5,
               n_ts=5,
               n_core=4,
               **options):
        if begin_time is None:
            begin_time = self.dataset.begin_time
        if end_time is None:
            end_time = self.dataset.end_time
        # endtime should be the last available time from dataset
        self.time = np.linspace(begin_time, end_time, n_ts)
        self.options.update(options)

        if n_core > 1:
            pool = mp.Pool(n_core)
            pool_results = []
            for i in range(n_iter):
                pool_results.append(
                    pool.apply_async(self.run,
                                     kwds={
                                         "begin_time": begin_time,
                                         "end_time": end_time,
                                         "n_iter": 1,
                                         "n_ts": n_ts,
                                         "mp": {
                                             "mp_i_iter": i,
                                             "n_core": n_core
                                         }
                                     }))
            pool.close()
            pool.join()
            self.list_of_parameter_trajectories = []
            self.list_of_state_trajectories = []
            for res_obj in pool_results:
                ptraj, straj = res_obj.get()
                self.list_of_parameter_trajectories.append(ptraj[0])
                self.list_of_state_trajectories.append(straj[0])

        else:
            self.run(begin_time=begin_time,
                     end_time=end_time,
                     n_iter=n_iter,
                     n_ts=n_ts)
        # convert the lists into xarray
        self.parameter_trajectories = xr.DataArray(
            data=np.array(self.list_of_parameter_trajectories),
            coords=[("iter", list(range(n_iter))), ("time", self.time),
                    ("param", list(self.parameter_names))],
            name="parameter trajectories")
        self.state_trajectories = xr.DataArray(
            data=np.array(self.list_of_state_trajectories),
            coords=[("iter", list(range(n_iter))), ("time", self.time),
                    ("state", list(self.dataset.get_state_names()))],
            name="state trajectories")
        self.parameter_trajectories.attrs.update(self.options)
        self.state_trajectories.attrs.update(self.options)
        return self.parameter_trajectories, self.state_trajectories

    def run(self, begin_time=None, end_time=None, n_iter=5, n_ts=5, **kw):
        self.list_of_parameter_trajectories = list()
        self.list_of_state_trajectories = list()

        for i_iter in range(n_iter):
            if self.options['verbose'] >= self.ITER:
                if "mp" in kw:
                    print(
                        f"iteration: {kw['mp']['mp_i_iter']} ({ mp.current_process().name })"
                    )
                else:
                    print(f"iteration: {i_iter}")
            self.parameter_trajectory = np.zeros(
                (n_ts, len(self.parameter_names)))
            self.state_trajectory = np.zeros((n_ts, len(self.dataset)))

            data = self.dataset.interpolate(n_ts=n_ts)
            # ! 👇 is probably wrong (*initial value problem*)
            # params = self.find_init_guesses()
            if self.options['init_method'] is None:
                self.state_trajectory[0, :] = data[:, 0, 0]
                self.parameter_trajectory[0, :] = self.parameters['init']
            elif self.options['init_method'] == "interfacefocus":
                # putting state initialization here for the possible return
                # of states in natal_init
                self.state_trajectory[0, :] = data[:, 0, 0]
                self.parameter_trajectory[0, :] = self.natal_init(i_iter, data)
            else:
                raise Exception("Unknown initialization method.")

            for i_ts in range(1, n_ts):
                if (i_ts %
                        10) == 0 and self.options['verbose'] >= self.TIMESTEP:
                    print(f"time step: {i_ts}")

                (
                    self.parameter_trajectory[i_ts, :],
                    self.state_trajectory[i_ts, :],
                ) = self.fit_timestep(
                    initial_guess=self.parameter_trajectory[i_ts - 1, :],
                    begin_states=self.state_trajectory[i_ts - 1, :],
                    interp_data=data[:, i_ts, :],
                    i_iter=i_iter,
                    i_ts=i_ts,
                )

            self.list_of_parameter_trajectories.append(
                self.parameter_trajectory)
            self.list_of_state_trajectories.append(self.state_trajectory)
        return self.list_of_parameter_trajectories, self.list_of_state_trajectories

    def natal_init(self, i_iter, data):
        """ Natal van Riel's parameter initialization
        Basically, it tries multiple times to fit the parameters to the data
        return: arrays of initial parameters
        """
        sse = np.inf
        i_ts = 0
        sseThres = self.options["sseThres"]

        while sse > sseThres:
            params = self.parameters['init'] * (10**(
                2 * np.random.random_sample(size=(len(self.parameter_names), ))
                - 1))
            lsq_res = least_squares(
                self.objective_function,
                params,
                bounds=(self.parameters["lb"], self.parameters["ub"]),
                method=self.options["method"],
                kwargs={
                    "begin_states": data[:, i_ts, 0],
                    "interp_data": data[:, i_ts, :],
                    "time_span": [-1000, 0],
                    "R": None,
                    "i_iter": i_iter,
                    "i_ts": i_ts
                },
            )
            sse = lsq_res.cost
        return lsq_res.x

    def steady_states(self):
        """ the data interpolation is assumed to be in steady states"""
        pass

    def fit_timestep(self,
                     initial_guess=None,
                     begin_states=None,
                     interp_data=None,
                     i_iter=None,
                     i_ts=None):
        """call least_squares
        access Optimizer options via `i_iter` and `i_ts` and `self`
        """
        time_span = self.time[i_ts - 1:i_ts + 1]
        bounds = (self.parameters["lb"], self.parameters["ub"])

        lsq_result = least_squares(
            self.objective_function,
            initial_guess,
            bounds=bounds,
            method=self.options["method"],
            kwargs={
                "begin_states": begin_states,
                "interp_data": interp_data,
                "time_span": time_span,
                "R": self.options['regularization'],
                "i_iter": i_iter,
                "i_ts": i_ts
            },
        )

        # compute the states using new parameters from lsq (a bit tedious yeah...)
        # but adding `begin_states` and `lsq_result.fun` can be confusing to users
        new_state_traj = self.model.compute_states(
            lsq_result.x,
            time_span,
            begin_states,
            new_param_names=self.parameter_names)
        return (lsq_result.x, new_state_traj[:, -1])

    def objective_function(self,
                           params,
                           begin_states=None,
                           interp_data=None,
                           time_span=None,
                           R=default_regularization,
                           i_iter=None,
                           i_ts=None,
                           **kw):
        """ Objective function

            The function minimized by least squares method. For ADAPT, an objective
            should:

            1. compute the states at the end of the `timespan`, using the give
                `parameters`, `begin_states`.
            2. choose those `end_states` and `interp_states` which are "observable"
            3. calculate and choose observable fluxes
            4. calculate residual
            5. calculate regularization term
            6. concatenate errors and regularizations

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
            R
                callable, regularization function. ❓: is it necessary to pass R
                as an argument?
            parameter_names
                list of strings, names of the parameters to optimize, other parameters
                are set as fixed (constant) parameters

            Return
            ------
            np.ndarray: shape(len(states) + len(parameter penalty))
        """

        end_states = self.model.compute_states(
            new_params=params,
            time_points=time_span,
            x0=begin_states,
            new_param_names=self.parameter_names)
        end_states = end_states[:, -1]
        # observable: choose those observable to compare with the data
        end_states = end_states[self.model.states["observable"]]

        # data doesn't contain unobservable states/fluxes, no mask
        interp_states = interp_data[:, 0]
        interp_stds = interp_data[:, 1]

        # equation 3.4 [ADAPT 2013]
        errors = (end_states - interp_states) / interp_stds
        """ need more work on the errors of flux. toy and trehalose don't have
        flux data, so I won't waste my time on it.
        New comer's job:
          - add flux in basemodel
          - add flux in dataset
          - add flux errors here """
        if R is not None:
            reg_term = self.options['lamda'] * R(
                params=params,
                parameter_trajectory=self.parameter_trajectory,
                state_trajectory=self.state_trajectory,
                i_iter=i_iter,
                i_ts=i_ts,
                time_span=time_span)
        else:
            reg_term = np.zeros((len(params)))  # SB
        residual = np.concatenate([errors, reg_term])
        return residual


def optimize(model,
             dataset,
             *params,
             n_iter=10,
             n_tstep=100,
             n_core=4,
             **options):
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
    optim.run_mp(n_iter=n_iter, n_ts=n_tstep, n_core=n_core, **options)
    return optim.parameter_trajectories, optim.state_trajectories, optim.time


def ADAPT():
    """ The ADAPT function
    add a threshold to elimimate the trajectories that don't fit enough.
    In other word, the parameters have to be random at the beginning.
    """
    # TODO adapt routine in interface focus code base, threshold


def steady_states(model, s):
    pass


if __name__ == "__main__":
    from pyADAPT.examples.lotka import LotkaVolterra
    from pyADAPT.examples.toy import ToyModel
    from pyADAPT.dataset import DataSet

    model = ToyModel()
    data = DataSet(
        raw_data_path="data/toyModel/toyData.mat",
        data_specs_path="data/toyModel/toyData.yaml",
    )
    ptraj, straj, time = optimize(model,
                                  data,
                                  "k1",
                                  n_iter=4,
                                  n_tstep=50,
                                  init_method=None,
                                  verbose=Optimizer.TIMESTEP)
