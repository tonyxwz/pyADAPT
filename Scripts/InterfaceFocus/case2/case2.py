# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:11:50 2020

@author: tonyx

case 2:Estimate lumped parameter k1 for each dataset
separately using a Monte Carlo approach.
Transform data provided as a mean value with a standard deviation into a
Gaussian distribution, from which samples are generated.
pEst is a matrix with k1 at the 5 different stages in the rows
"""

import numpy as np
import matplotlib.pyplot as plt

from ode import ode
from fluxes import fluxes


    