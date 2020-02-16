""" optimize module is a *procedure* version of ADAPT

_ADAPT I am considering renaming ADAPT to Optimizer in order to be distant away from
this lame paper.
"""
import datetime
import multiprocessing as mp
import os
import sys

import numpy as np
import pandas as pd
from lmfit import Parameter, minimize
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import least_squares, leastsq

from .dataset import DataSet
from .model import Model

#%% define ADAPT outline function


def optimize():
    """ the main optimization (ADAPT) procedure 
    
    """
    # 1. randomize
    # 2.
