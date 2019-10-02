# -*- coding: utf-8 -*-
import logging

import numpy as np
from lmfit import Parameters, minimize
from scipy.integrate import OdeSolver, odeint, solve_ivp
from scipy.optimize import least_squares

from pyADAPT import DataSet, Model


TOLERANCE = 100

class ADAPTOptions(object):
    """ helper """
    pass

class ADAPTapp(object):
    """ implementation of ADAPT algorithm 
    Encapsulates all the methods required to perform ADAPT simulation
    """
    def __init__(self, model: Model, dataset: DataSet):
        """ Design paradigm: the model and dataset should not store any
        runtime data, ADAPTapp instance handles everything else
        """
        # *1 logger
        self.logger = logging.getLogger(self.__str__())
        self.logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARNING)
        fh = logging.FileHandler("adapt.log")
        fh.setLevel(logging.INFO)

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.addHandler(fh)

        # *2 simulation options
        self.options: dict = dict()
        self.options['n_ts']=  100  # number of time steps
        self.options['n_iter'] = 20  # number of iteration
        self.options['smin'] = -1  # scale min
        self.options['smax'] = 1  # scale max
        self.options['seed'] = 0  # to seed the random generator
        self.options['ss_time'] = 100  # steady-state duration
        

        # *3 model components
        self.model = model
        self.dataset = dataset
        # TODO check the names consistency in model and dataset

        # *4 simulation results
        self.current_timestep = 0
        self.trajectories:list = list()  # index being the iteration number
        # trajectories[0]['s1'][0]

    def run(self, n_iter=None):
        np.random.seed(self.options['seed'])
        if n_iter is None:
            n_iter = self.options['n_iter']
        for i in range(n_iter):
            self.trajectories.append(dict())
            # self.logger.info('Start iteration %d' % i)
            # 1. randomize data (generate splines)
            # d is the randomized data, which is a dict {s1, s2} of dict {value, std}
            d = self.randomize_data()
            # 2. randomize parameter set (pagina 32 op thesis)
            # p: initial values of the parameters, dict {s1, s2}
            p = self.randomize_parameters()
            # 3. optimize the penalty function
            self.steady_state_sim(duration=100)
            for j in range(len(np.arange(10))):  # TODO change to time step no.
                # self.logger.info('Fit time step %d' % j)
                self.fit_timestep(j)

    def preintervention_sim(self):
        pass

    def steady_state_sim(self, duration=100):
        """ steady state simulation before intervention (t0)
        `p` at the end this simulation is kept as the initial guesses
        """
        pass

    def fit_timestep(self, tstep, phase="fit"):
        """ simulate one timestep.
        1. compute states
        2. compute fluxes
        """
        # self.logger.error('not implemented')
        if phase == "fit":
            pass
        elif phase == "ss": # steady state
            pass
        else:
            raise Exception(f"Unknown option: {phase}")

    def optimize_parameters(self):
        """ 1. find which parameters are fittable
        2. optimize with scipy.optimize.least_squares. The differential
        equations need be solved in the objective function passed into the
        optimizer.
        TODO: define the objective function.
        """
        pass
        # scipy.optimize.least_squares

    @staticmethod
    def objective_func(p, x0, tspan, tstep):
        """ Objective function that gives the error between model output and
        the penalty for parameter changes.

        Parameters
        ----------
        `p` : list
            parameters to evaluate the error with

        `x0` : numpy.ndarray
            initial states as the initial value of the ode

        `tspan` : numpy.ndarray of size (2,)
            time span to integrate on

        `tstep` : int
            the ID of the current time step
        """
        pass
        # scipy.integrate.solve_ivp

    def randomize_data(self):
        d = self.dataset.generate_splines(self.options['n_ts'])
        return d

    def randomize_parameters(self):
        """randomize the parameters at the beginning of a simulation
        using the formula in van Beek's thesis and matlab function in:
        *AMF.Model.randomizeParameters*
        """
        p = self.model.randomize_params(self.options['smin'], self.options['smax'])
        for k, v in p.items():
            self.trajectories[-1][k] = [v]  # always handling the newest traj
        return p

    def regularization_term(self, p):
        # theta: new parameters
        # Xi^2(theta[n]) = sum((theta[n][i] - theta[n-1][i]) /
        #                       delta_t * 1 / theta[0][i])
        # from i=1 tot N_p (number of parameters)
        p_init = self.trajectories[0]
        p_previous = self.trajectories[-1]
        penalty = np.sum(((p - p_previous) / self.delta_t * 1 / p_init)**2)
        return penalty

    def set_options(self, **kwargs):
        for k, v in kwargs.items():
            if k in self.options:
                self.options[k] = v
            else:
                self.logger.error(f"Unknown option: {k}")
