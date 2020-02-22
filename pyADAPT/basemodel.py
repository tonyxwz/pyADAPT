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

TODO
====
1. warning about the unused names
2. clever way to add ofi
"""
# import networkx as nx
from collections import OrderedDict
from abc import abstractmethod, ABCMeta

import numpy as np
from scipy.integrate import solve_ivp
from lmfit import Parameters, minimize, Parameter

from pyADAPT.io import read_data_specs


class BaseModel(metaclass=ABCMeta):
    """Abstract class for constructing ADAPT models

    One should inherit from this class and extend the constructor method and the
    `ode` method to specify the structure of the model

    1. Define the model parameters, constants and states in the constructor method.
    2. Define the ODE in `odefunc`, please read `lmfit` documentation:
        https://lmfit.github.io/lmfit-py/parameters.html.
    3. Define model input function.
    4. Define reactions function.
    """
    def __new__(cls, *args, **kwargs):
        # to use `super` in `__new__` method: issubclass(cls, Model) is True
        # super() is (object or type) args and kwargs should be handle in this method
        instance = super().__new__(cls)  # , *args, **kwargs)  # model instance
        instance.name = "Base Model"
        instance.notes = " ".join([
            "This model should not be instantiated for it only",
            "serves as the base class of other models, please refer to the",
            "docstring about how to extend this class and define your own",
            "model.",
        ])
        instance.specs = {}  # TODO
        # time, named as convention from matlab version
        # TODO change predictor as a component of ADAPT
        # with a default value, if cannot be deduced from the Dataset
        instance.predictor = []

        # ? why OrderedDict
        # ! Ans: OrderedDict guarrantees that the order of the dictionary and the order
        # ! of the array e.g. observables recorded are the same. (odict.values())
        instance.constants = OrderedDict()
        instance.parameters = Parameters()
        instance._init_parameters = OrderedDict()
        # instance.reactions = list()
        instance.states = OrderedDict()
        instance.init_vary = []
        instance.observables = OrderedDict()
        instance.state_musk = list()  # observable state
        instance.flux_musk = list()  # observabel flux

        instance.variable_names = set()
        # instance.i_tstep = 0
        # instance.i_iter = 0
        return instance

    def __init__(self):
        """REMARK: The differences between a parameter that will not be fitted
        " and a constant in the model:
        "
        " the parameters are still possible to be optimized by minimizer but the
        " constant are not. This provides the flexibility to select which
        " parameters to fit.
        """
        self.name: str
        self.notes: str
        self.spec: dict
        self.predictor: list
        self.constants: OrderedDict
        self.parameters: Parameters
        self._init_parameters: OrderedDict
        self.states: OrderedDict
        self.init_vary: list
        self.observables: OrderedDict
        self.variable_names: set
        self.state_musk: list
        self.flux_musk: list
        self.vary_flags = [p.vary for p in self.parameters.values()]

    # @staticmethod
    @abstractmethod
    def fluxes(self, t, x, p):
        raise NotImplementedError

    # @abstractmethod
    def inputs(self, t):
        raise NotImplementedError

    # @staticmethod
    @abstractmethod
    def odefunc(self, t, x, p):
        """ t: time, x: state, p: parameters, u: constant """
        raise NotImplementedError

    def regfunc(self, p, p_prev, p_init, delta_t):
        """ regfunc is not a abstract method because it provides a
        default regularization.
        """
        pass

    def compute_states(self,
                       t_span,
                       x0,
                       p=None,
                       rtol=1e-7,
                       atol=1e-7,
                       t_eval=None):
        """ no `uvec`: given time `t`, input should be determined as well.
        " From here, the odefunc is integrated to get the state variables.
        " Then state variables will be used to evaluate the reactions, and the
        " reactions will be compared against the experimental data to minimize
        " parameters.
        """
        if p is None:
            p = self.parameters
        if t_eval is None:
            t_eval = [t_span[-1]]
        # u = self.inputs(t_span[0])  # input
        if type(x0) is OrderedDict or type(x0) is dict:
            x0 = np.array(list(x0.values()))

        sol = solve_ivp(
            lambda t, x: self.odefunc(t, x, p),
            t_span,
            x0,
            t_eval=t_eval,
            rtol=rtol,
            atol=atol,
        )

        # solve_ivp always return an array of shape(len(x0), len(t_eval))
        # so I have to squeeze the result if I want only the last row
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        # y is in the same order as x0
        return sol.y

    def compute_reactions(self, t, x, p):
        # ? necessary to be in consistent with the matlab api?
        return self.fluxes(t, x, p)

    def add_parameter(
            self,
            name="",
            value=None,
            vary=True,
            lb=-np.inf,
            ub=np.inf,
            expr=None,
            brute_step=None,
    ) -> None:
        """Read the docs here:
        https://lmfit.github.io/lmfit-py/parameters.html#the-parameter-class
        """
        self.add_name(name)
        self.parameters.add(
            name=name,
            value=value,
            vary=vary,
            min=lb,
            max=ub,
            expr=expr,
            brute_step=brute_step,
        )
        self._init_parameters[name] = value

    def add_constant(self, name="", value=None) -> None:
        self.add_name(name)
        self.constants[name] = value

    def add_predictor(self, name="", value=[]) -> None:
        self.add_name(name)
        self.predictor = value
        self.predictor_name = name

    def add_state(self, name="", init=None, observable=True) -> None:
        self.add_name(name)
        self.states[name] = init
        self.state_musk.append(observable)

    def randomize_params(self, smin, smax, params=None):
        """using the formula in van Beek's thesis and matlab function in:
        `AMF.Model.randomizeParameters`

        smin, smax
        """
        p = self.parameters.copy()
        # if params is None:
        #     params = self.init_vary
        for k, v in p.items():  # TODO maybe change v.vary to initvary flags
            if v.vary:  # values of dict observables is boolean
                v.value = p[k] * np.power(10, (
                    (smax - smin) * np.random.rand() + smin))
        return p

    def add_name(self, name: str):
        # TODO: add reference to the variable
        if self.check_name(name):
            self.variable_names.add(name)
        else:
            raise Exception("%s is already define in the model" % name)

    def check_name(self, name):
        return not name in self.variable_names


if __name__ == "__main__":
    pass
