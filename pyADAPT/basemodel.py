"""No configuration files for models anymore except for a python script that extends this
abstract class `Model`.

Steps
=====
1: toymodel without considering the model input
        with only constant parameters               *
3: toymodel validation
4: Tiemann's model, fully implemented
5: Clamp model in Pascal van Beek's work
6: Discuss the Model with David and Natal
"""
from abc import ABCMeta, abstractmethod
from collections import OrderedDict

# import networkx as nx
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from pyADAPT.io import read_data_specs


class BaseModel(metaclass=ABCMeta):
    """Abstract class for constructing ADAPT models

    One should inherit from this class and extend the constructor method and the
    `ode` method to specify the structure of the model
    """
    def __new__(cls, *args, **kwargs):
        instance = super().__new__(cls)
        instance.name = "Base Model"
        instance.notes = """Base Model, for extending only"""
        instance.predictor_name = "time"
        instance.map = dict()
        instance.fluxes = list()
        instance._parameters = list()
        instance._states = list()
        instance._constants = list()

        return instance

    def __init__(self):
        self.name: str
        self.map: dict
        self.fluxes: list
        self._parameters: list
        self._states: list
        self._constants: list

        # not checking symbol names
        self.parameters = pd.DataFrame(
            self._parameters,
            columns=["name", "value", "vary", "lb", "ub", "init"],
            index=[p[0] for p in self._parameters],
        )
        # only init values of states matters in the model
        self.states = pd.DataFrame(
            self._states,  # value in Model.states are not used in computation
            columns=["name", "value", "observable", "init"],
            index=[s[0] for s in self._states],
        )
        self.constants = pd.DataFrame(self._constants,
                                      columns=['name', 'value'],
                                      index=[c[0] for c in self._constants])
        del self._parameters
        del self._states
        del self._constants

    def add_parameter(self,
                      name="",
                      value=None,
                      vary=False,
                      lb=-np.inf,
                      ub=np.inf,
                      **kw):
        # name, value, vary?, lb, ub, init
        self.map[name] = len(self._parameters)
        self._parameters.append([name, value, vary, lb, ub, value])

    def add_parameters(self, ll):
        for r in ll:
            if type(r) is list:
                self.add_parameter(*r)
            elif type(r) is dict:
                self.add_parameter(**r)

    def add_state(self, name='', value=None, observable=False):
        # name, value, observable, init
        self.map[name] = len(self._states)
        self._states.append([name, value, observable, value])

    def add_states(self, ll):
        for r in ll:
            if type(r) is list:
                self.add_state(*r)
            elif type(r) is dict:
                self.add_state(**r)

    def add_constant(self, name='', value=None):
        self._constants.append([name, value])

    def add_constants(self, ll):
        for r in ll:
            if type(r) is list:
                self.add_constant(*r)
            elif type(r) is dict:
                self.add_constant(**r)

    @abstractmethod
    def ode_func(self, t, x, p):
        """ t: time
        x: state at t, ndarray
        p: parameters, DataFrame
        return: dx/dt at t

        TODO
        ----
        further simplify the process for indexing states x and parameters p from
        the arguments.
        """
        pass

    def jacobian_func(self, t, x, p):
        # This might be useful if ODE solvers require it.
        raise NotImplementedError

    def fluxes_func(self, t, x, p):
        raise NotImplementedError

    def inputs_func(self, t):
        raise NotImplementedError

    def __repr__(self):
        return " ".join([super().__repr__(), self.name])

    def compute_states(
        self,
        new_params=[],  # parameters that need to be optimized
        time_points=None,  # the time span of the computation
        x0=None,  # the states at the first time point
        new_param_names=[],  # the parameter's names, in the same order
        method="RK45",  # odesolver, only RK45/RK23/DOP853 since no jacobians
        rtol=1e-3,  # relative tolerance
        atol=1e-6  # absolute tolerance
    ):
        if new_param_names:
            # https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#why-does-assignment-fail-when-using-chained-indexing
            # see the `def do_something(df)` example
            self.parameters.loc[new_param_names, "value"] = new_params
            # parameters.loc[new_param_names] = new_params

        t_span = [time_points[0], time_points[-1]]

        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        sol = solve_ivp(
            lambda t, x: self.ode_func(t, x, self.parameters['value']),
            t_span,
            x0,
            method=method,
            t_eval=time_points,
            rtol=rtol,
            atol=atol,
        )

        return sol.y

    # move this method to Optimizer
    def randomize_params(self, smin, smax, params=None):
        """using the formula in van Beek's thesis and matlab function in:
        `AMF.Model.randomizeParameters`

        smin, smax
        """
        p = self.parameters.copy()
        # if params is None:
        #     params = self.init_vary
        for k, v in p.items():  # maybe change v.vary to initvary flags
            if v.vary:  # values of dict observables is boolean
                v.value = p[k] * np.power(10, (
                    (smax - smin) * np.random.rand() + smin))
        return p

    def sync_from_optimizer(self, optimizer):
        """ called at the end of Optimizer.run
        """
        for i in range(len(optimizer)):
            self.parameters.iloc[i]["value"] = optimizer[i]

    def reset(self):
        """ reset to initial conditions """
        self.parameters['value'] = self.parameters['init']
        self.states['value'] = self.states['init']

    def draw(self):
        """ Draw the model using networkx """


if __name__ == "__main__":
    pass
