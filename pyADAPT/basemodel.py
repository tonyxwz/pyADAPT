# -*- coding: utf-8 -*-
"""No configuration files for models anymore except for a python script that extends this
abstract class `Model`.

Steps
=====
1: toymodel without considering the model input
        with only constant parameters ... ok
2: toymodel validation ... ok
3: Discuss the Smallbone Model with David and Natal
"""
from abc import ABCMeta, abstractmethod
from collections import OrderedDict

# import networkx as nx  # visualize the Model as a network
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from pyADAPT.io import read_data_specs


class ParameterTable(pd.DataFrame):
    def __init__(self):
        super().__init__()

    def __iadd__(self, other):
        pass

    def __add__(self, other):
        pass


class BaseModel(metaclass=ABCMeta):
    """Abstract class for constructing ADAPT models

    One should inherit from this class and extend the constructor method and the
    `ode` method to specify the structure of the model
    """

    def __new__(cls, *args, **kwargs):
        instance = super().__new__(cls)
        instance.name = "Base Model"
        instance.notes = """Base Model, for extending only"""

        instance._parameters = list()
        instance.map = dict()

        # instance.mat = np.ndarray()  # stoichiometry matrix
        return instance

    def __init__(
        self, state_order=[], initial_states=[], flux_order=[], input_order=[]
    ):
        self.name: str
        self.notes: str
        self.map: dict
        self.mat: np.ndarray

        self.state_order = state_order

        self.initial_states = np.array(initial_states)
        self.flux_order = flux_order
        self.input_order = input_order

        self._parameters: list
        self.parameters = pd.DataFrame(
            self._parameters,
            columns=["name", "value", "vary", "lb", "ub", "init"],
            index=[p[0] for p in self._parameters],
        )
        del self._parameters

        # ? add map for parameters
        for i, s in enumerate(self.state_order):
            self.map[s] = i
            setattr(self, s, i)
        for i, v in enumerate(self.flux_order):
            self.map[v] = i
            setattr(self, v, i)
        for i, ip in enumerate(self.input_order):
            self.map[ip] = i
            setattr(self, ip, i)

    def __getitem__(self, x):
        return self.map[x]

    def add_parameter(self, name, value, vary, lb=-np.inf, ub=np.inf, **kw):
        # name, value, vary?, lb, ub, init
        self._parameters.append([name, value, vary, lb, ub, value])
        self.map[name] = len(self._parameters)

    def add_parameters(self, l: list):
        for p in l:
            self.add_parameter(**p)

    @property
    def has_flux(self):
        return len(self.flux_order) > 0

    @abstractmethod
    def state_ode(self, t, x, p):
        """ t: time, x: state at t, p: parameters, u: constant
        return: dx/dt at t
        """
        pass

    def fluxes(self, t, x, p):
        raise NotImplementedError

    def inputs(self, t):
        raise NotImplementedError

    def __repr__(self):
        return " ".join([super().__repr__(), self.name])

    def compute_states(
        self,
        time_points=None,  # the time span of the computation
        x0=None,  # the states at the first time point
        new_params=[],  # parameters that need to be optimized
        new_param_names=[],  # the parameter's names, in the same order
        odesolver="RK45",  # odesolver
        rtol=1e-3,  # relative tolerance TODO
        atol=1e-6,  # absolute tolerance
        **solver_options
    ):
        # integrate state_ode to get state values
        if new_param_names:
            self.parameters.loc[new_param_names, "value"] = new_params

        # t_span = np.linspace(time_points[0], time_points[-1], 1000)
        t_span = [time_points[0], time_points[-1]]
        sol = solve_ivp(
            lambda t, x: self.state_ode(t, x, self.parameters["value"]),
            t_span,
            x0,
            method=odesolver,
            t_eval=time_points,
            rtol=rtol,
            atol=atol,
            **solver_options
        )
        return sol.y

    # move this method to Optimizer
    def randomize_params(self, smin, smax, params=None):
        """using the formula in van Beek's thesis and matlab function in:
        `AMF.Model.randomizeParameters`

        smin, smax
        """
        # TODO remove
        p = self.parameters.copy()
        # if params is None:
        #     params = self.init_vary
        for k, v in p.items():
            if v.vary:  # values of dict observables is boolean
                v.value = p[k] * np.power(10, ((smax - smin) * np.random.rand() + smin))
        return p

    def sync_from_optimizer(self, optimizer):
        """ called at the end of Optimizer.run
        """
        for i in range(len(optimizer)):
            self.parameters.iloc[i]["value"] = optimizer[i]

    def reset(self):
        """ reset to initial conditions """
        self.parameters["value"] = self.parameters["init"]


if __name__ == "__main__":
    pass
