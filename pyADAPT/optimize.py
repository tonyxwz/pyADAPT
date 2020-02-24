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

[^1]: I am considering renaming ADAPT to Optimizer in order to be distant away
from this lame paper.
"""

import datetime
import multiprocessing as mp  # consider dask
import os
import sys

import numpy as np
import pandas as pd
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import least_squares, leastsq

# ? tentative, as model and data set should be pass as arguments to the function in this module, these will only be used as type hint
from pyADAPT.dataset import DataSet
from pyADAPT.basemodel import BaseModel


class Optimizer(object):
    """optimizes an ADAPT model"""

    def __init__(self, model, dataset):
        self.model = model
        self.dataset = dataset
        self.lamda = 1  # weight of the regularization

    def run(self):
        pass

    def report(self):
        pass

    def find_init_guesses(self):
        """ Find the initial guess of the parameters at the start of each iteration.
        Problem statement: the data need to be randomized at time 0. we would like to find a set of parameters that will lead to a steady states of the randomized states at t0.

        # TODO solve problem above
        """


def optimize(model, dataset):
    """ the main optimization (ADAPT) procedure
    """
    # 1. randomize
    # 2.


# TODO move to Optimizer
def objective_function(params, model: BaseModel, x_begin, d_end, time_span, R, L):
    """Objective function

    The function that is called by the minimizers. In the case of ADAPT, the
    objective function should first simulate the model from the begin of
    `time_span` till the end, using `x_begin` as the initial state (ivp). After
    the simulation, the ode solver should return an array of states. We only
    want the last column (snapshot) to be compared with the experimental data
    `d_end`.

    Parameters
    ----------
    params
        the parameters θ of the model
    model
        pyADAPT.model.Model instance
    x_begin
        the initial states
    d_end
        the experimental data at the end of the time span, containing the
        states and corresponding standard deviations.
    time_span
        the time span to solve the ODE
    R
        callable, regularization function
    L
        lambda, the weight of the regularization term
    """
    states = model.compute_states(time_span, x_begin, params)
    x_end = states[:, -1]
    # select those observable/measurable states to compare with the data
    ox_end = x_end[model.state_musk]
    # TODO use pandas for array musking
    s = d_end[model.state_musk, 0]
    v = d_end[model.state_musk, 1]

    # ? I think it only makes sense in fake data that the data interpolation
    # ? can contain unmeasurable data. otherwise, the following should be used.
    # ? assertion is needed to check ox_end, s and v have the same dimension.
    # s = d_end[0]
    # v = d_end[1]

    # ? equation 3.4 [ADAPT 2013]
    errors = (ox_end - s) / v
    reg_term = L * R()  # TODO parameters of `R`
    residual = np.concatenate([errors, reg_term])
    return residual


def reg_example(p=None,):
    # TODO find literatures about regularization, to replace eq3.5
    pass


def steady_states(model, s):
    pass


if __name__ == "__main__":
    pass
