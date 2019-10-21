# -*- coding: utf-8 -*-
import logging
from multiprocessing import Pool, current_process

import numpy as np
from lmfit import Parameters, minimize
from scipy.integrate import OdeSolver, odeint, solve_ivp
from scipy.optimize import least_squares
from cached_property import cached_property

from pyADAPT import DataSet, Model



class ADAPTOptions(object):
    """ helper """
    pass

class ADAPT(object):
    """ implementation of ADAPT algorithm
    Encapsulates all the methods required to perform ADAPT simulation
    """
    def __init__(self, model: Model, dataset: DataSet):
        """ Design paradigm: the model and dataset should not store any
        runtime data, ADAPT instance handles everything else

        In `ADAPT` class, the identifier of a name should be key of the dictionary
        But in class `DataSet` and class `Model`, musks and flags should be used.
        """
        # *1 logger
        self.logger = logging.getLogger('ADAPT')
        self.logger.setLevel(logging.WARNING)
        ch = logging.StreamHandler()
        # ch.setLevel(logging.WARNING)
        fh = logging.FileHandler("adapt.log")
        # fh.setLevel(logging.INFO)

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.addHandler(fh)

        # *2 default simulation options
        self.options: dict = dict()
        self.options['n_ts']=  100  # number of time steps
        self.options['n_iter'] = 20  # number of iteration
        self.options['smin'] = -1  # scale min
        self.options['smax'] = 1  # scale max
        self.options['seed'] = 0  # to seed the random generator
        self.options['ss_time'] = 100  # steady-state duration
        self.options['lambda'] = 0.1  # regularization weight
        self.options['rtol'] = 1e-7
        self.options['atol'] = 1e-7
        self.options['log_level'] = logging.WARNING


        # *3 model components
        self.model = model
        self.dataset = dataset
        # TODO check the names consistency in model and dataset

        # *4 preparation
        # *5 simulation results
        self.i_iter = 0
        self.i_tstep = 0
        self.time_points:list = []
        self.current_timestep = 0
        self.states = np.zeros((1,))
        self.trajectories = np.zeros((1,))  # index being the iteration number
        self.min_history:list = []


    @cached_property
    def delta_t(self):
        return (self.model.predictor[-1] - self.model.predictor[0]) / self.options['n_ts']

    def set_options(self, **kwargs):
        for k, v in kwargs.items():
            if k in self.options:
                self.options[k] = v
            else:
                self.logger.error(f"Unknown option: {k}")
                raise Exception(f"Unknown option: {k}")

    def get_tspan(self, ts):
        # time range of timestep ts
        if ts == 0:
            tspan = [ self.model.predictor[0]-self.options['ss_time'],
                    self.model.predictor[0] ]
        else:
            tspan = self.time_points[ts-1: ts+1]
        return tspan

    def add_to_trajectories(self, p=None, i_iter=None, i_tstep=None):
        # append the parameters to the latest trajectory (current iteration)
        if p is None:
            p = np.array(list(self.model.parameters.values()))
        if i_iter is None:
            i_iter = self.i_iter
        if i_tstep is None:
            i_tstep = self.i_tstep + 1
        self.trajectories[i_iter, :, i_tstep] = p

    def randomize_data(self):
        d = self.dataset.interpolate(self.options['n_ts'])
        return d

    def randomize_init_parameters(self, append=True):
        """randomize the parameters at the beginning of a simulation
        """
        self.model.randomize_params(self.options['smin'], self.options['smax'])
        # always at tstep=0
        if append:
            self.add_to_trajectories(i_tstep=0)

    def run(self, n_iter=None, n_core=0):
        np.random.seed(self.options['seed'])

        if n_iter is None:
            n_iter = self.options['n_iter']
        else:
            self.options['n_iter'] = n_iter

        self.logger.setLevel(self.options['log_level'])
        self.time_points = np.linspace(self.model.predictor[0], self.model.predictor[1], self.options['n_ts'])

        self.trajectories = np.zeros((n_iter, len(self.model.parameters), self.options['n_ts']+1))
        self.states = np.zeros((n_iter, len(self.model.states), self.options['n_ts']+1))

        if n_core > 1:
            pool = Pool(n_core)
            pool_results = []
            # for i_iter in range(n_iter):
            #     r_ = pool.apply_async(self.parallel_kernel, 
            #             args=(i_iter,

            #             )
            #         )
            #     pool_results.append(r_)
            pool.close()
            pool.join()
            # for r_ in pool_results:
            #     mtx = r_.get()
        else:
            for i_iter in range(n_iter):
                self.i_iter = i_iter
                self.logger.debug(f"iteration {i_iter}")
                # 1. randomize data (generate splines)
                # d [n step]
                data = self.randomize_data()
                # 2. randomize parameter (in place)
                self.randomize_init_parameters()
                self.states[i_iter, :, 0] = np.array(list(self.model.states.values()))
                for i_tstep in range(self.options['n_ts']):
                    # self.logger.debug(f"timestep {ts}")
                    self.i_tstep = i_tstep
                    d = data[:, i_tstep, :]  # select all the data at ts
                    # FIXME x0 should be the last moment of previous time step
                    # `self.model.states` should always be the initial states, `i_step` because init 
                    x0 = self.states[i_iter, :, i_tstep]
                    min_res = self.fit_timestep(x0, d, i_iter, i_tstep)
                    # FIXME min_history has no iteration dimension
                    self.min_history.append(min_res)
                    self.add_to_trajectories()

            #TODO return AdaptedModel

    def parallel_kernel(self, i_iter):
        self.logger.debug(f"{current_process()}: iteration {i_iter}")
        trajectory = np.zeros((len(self.model.parameters), self.options['n_ts']+1))
        states = np.zeros((len(self.model.states), self.options['n_ts']+1))
        data = self.randomize_data()
        self.randomize_init_parameters(append=False)
        trajectory[:, 0] = np.array( list(self.model.parameters.values()) )
        states[:, 0] = np.array( list(self.model.states.values()) )
        for i_tstep in range(self.options['n_ts']):
            # self.logger.debug(f"timestep {ts}")
            self.i_tstep = i_tstep
            d = data[:, i_tstep, :]  # select all the data at i_tstep
            x0 = self.states[i_iter, :, i_tstep]
            _min_res = self.fit_timestep(x0, d, i_iter, i_tstep)
            trajectory[:, i_tstep+1] = np.array( list(self.model.parameters.values()) )
            states[:, i_tstep+1] = np.array( list(self.model.states.values()) )

        return trajectory, states
        

    def fit_timestep(self, x0, d, i_iter, i_tstep):
        # x0 is the simulation result from previous time span
        # d is the interp data in current iteration
        # call objective function in minimizer here:
        t_span = self.get_tspan(i_tstep)

        minimization_result = minimize(self.objective_func, self.model.parameters,
                args=[x0, d, t_span, i_iter, i_tstep])
        # FIXME model should be static in this project, also i_iter, i_tstep
        self.model.parameters = minimization_result.params
        if not minimization_result.success:
            self.logger.warning(f'unsuccessful minimization, iter: {i_iter}, t_step: {i_tstep}')
        return minimization_result


    def objective_func(self, p, x0, d, t_span, i_iter, i_tstep):
        # 1. solve the ode using `p` and initial condition `x0`
        states = self.model.compute_states(t_span, x0, p,
                rtol=self.options['rtol'],
                atol=self.options['atol'],
                t_eval=[t_span[-1]])
        states = states[:, -1]  # select the last moment
        self.states[i_iter, :, i_tstep+1] = states
        f = self.model.compute_reactions(t_span[-1], states, p)
        # select observables
        observable_states = states[self.model.state_musk]
        observable_states_values = d[self.model.state_musk, 0]  # observable states from data (value)
        observable_states_stds = d[self.model.state_musk, 1]  # observable states from data (std)

        # FIXME fluxes are not computed
        # of = f[self.model.ofi]
        errors = (observable_states - observable_states_values) / observable_states_stds
        reg = self.tiemann_regularization(p, i_iter, i_tstep)
        errors = np.concatenate([errors, reg])
        # self.logger.debug(errors)
        return errors

    def tiemann_regularization(self, p, i_iter, i_tstep):
        """
        FIXME: regularization should be provided by the model
        theta: new parameters
        Xi^2(theta[n]) = sum((theta[n][i] - theta[n-1][i]) /
                              delta_t * 1 / theta[0][i])
        from i=1 tot N_p (number of parameters)
        """
        p_init = self.trajectories[i_iter, :, 0]
        p_previous = self.trajectories[i_iter, :, i_tstep]
        penalty = self.options['lambda'] * ((p - p_previous) / self.delta_t / p_init)**2
        return penalty

    def mkplot(self):
        # TODO aggregate several plotting methods in here
        pass