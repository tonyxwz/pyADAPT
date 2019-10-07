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
2. 
"""
# import networkx as nx
from collections import OrderedDict
from abc import abstractmethod, ABCMeta

import numpy as np
from scipy.integrate import solve_ivp
from lmfit import Parameters, minimize, Parameter

from .io import read_data_info


class ModelComponents(OrderedDict):
    # TODO fancy printing
    def __repr__(self):
        return super().__repr__()


class OverrideProtection(ABCMeta):
    """metaclass to protect private methods in `Model`"""
    def __new__(cls, name, bases, nmspc, **kwargs):
        protected_names = ['add_parameter', 'add_state',
                'add_constant', 'add_predictor']

        overriding = any([(x in nmspc) for x in protected_names])
        if bases in nmspc and overriding:
            raise Exception("Overriding protected components.")
        return type.__new__(cls, name, bases, nmspc)


def protect(*protected):
    """Returns a metaclass that protects all attributes given as strings"""
    class Protect(type):
        has_base = False
        def __new__(meta, name, bases, attrs):
            if meta.has_base:
                for attribute in attrs:
                    if attribute in protected:
                        raise AttributeError('Overriding of attribute "%s" not allowed.' % attribute)
            meta.has_base = True
            klass = super().__new__(meta, name, bases, attrs)
            return klass
    return Protect


class Model(metaclass=ABCMeta):
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
        instance = super().__new__(cls, *args, **kwargs)  # model instance
        instance.name = "Base Model"
        instance.description = " ".join([
            "This model should not be instantiated for it only",
            "serves as the base class of other models, please refer to the",
            "docstring about how to extend this class and define your own",
            "model."
        ])
        instance.specs = {}  # TODO
        # time, named as convention from matlab version
        instance.predictor = []

        # ? why OrderedDict
        # ! Ans: OrderedDict guarrantees that the order of the dictionary and the order
        # ! of the array e.g. observables recorded are the same. (odict.values())
        instance.constants = OrderedDict()
        instance.parameters = Parameters()
        # instance.reactions = list()
        instance.states = OrderedDict()
        instance.observables = OrderedDict()
        instance.trajectories = OrderedDict()
        
        instance.var_names = set()
        return instance

    def __init__(self):
        """REMARK: The differences between a parameter that will not be fitted 
        " and a constant in the model:
        "
        " the parameters are still possible to be optimized by minimizer but the
        " constant are not. This provides the flexibility to select which
        " parameters to fit.
        """
        self.n_tstep = 0  # the number of timesteps, initially 0
        self.n_iter = 0  # the number of iteration in the optimizer

    # @staticmethod
    @abstractmethod
    def reactions(self, t, x, p):
        raise NotImplementedError

    @abstractmethod
    def inputs(self, t):
        raise NotImplementedError

    # @staticmethod
    @abstractmethod
    def odefunc(self, t, x, p):
        """ t: time, x: state, p: parameters, u: constant """
        raise NotImplementedError

    def compute_states(self, t_span, x0, p=None,
            rtol=1e-7, atol=1e-7, t_eval=None):
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
            x0 = list(x0.values())

        sol = solve_ivp(lambda t, x: self.odefunc(t, x, p), t_span, x0,
                        t_eval=t_eval, rtol=rtol, atol=atol)

        # y is in the same order as x0
        return sol.y

    def compute_reactions(self, t, x, p):
        # ? necessary
        return self.reactions(t, x, p)

    def add_parameter(self, name="", value=None, vary=True,
            lb=-np.inf, ub=np.inf, expr=None, brute_step=None) -> None:
        """Read the docs here:
        https://lmfit.github.io/lmfit-py/parameters.html#the-parameter-class
        """
        self.add_name(name)
        self.parameters.add(name=name, value=value, vary=vary, min=lb, max=ub,
            expr=expr, brute_step=brute_step)

    def add_constant(self, name="", value=None) -> None:
        self.add_name(name)
        self.constants[name] = value

    def add_predictor(self, name='', value=[]) -> None:
        self.add_name(name)
        self.predictor = value
        self.predictor_name = name

    def add_state(self, name="", init=None, observable=True) -> None:
        self.add_name(name)
        self.states[name] = init
        self.observables[name] = observable

    def randomize_params(self, smin, smax):
        """smin, smax"""
        rp = OrderedDict()  # maybe dict is as good in python 3.7+
        # !!! Wrong !!!
        # TODO fix
        for k,v in self.observables.items():
            if v: # values of dict observables is boolean
                rp[k] = np.power(10, ((smax - smin) * np.random.rand() + smin))
        return rp

    def add_name(self, name: str):
        # TODO: add reference to the variable
        if self.check_name(name):
            self.var_names.add(name)
        else:
            raise Exception('%s is already define in the model' % name)

    def check_name(self, name):
        return not name in self.var_names

    def begin_extend(self):
        # do it in __new__
        pass

    def end_extend(self):
        # do it in __init__
        pass
