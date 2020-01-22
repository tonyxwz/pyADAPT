"""
case 2 and case 3 follow the same model structure. case 3 is a 
dynamic version of case 3.

F1 in the model of case 1 is simplified from

               Vmax   
    f1 =  ( --------- ) * u1(t) * S2
             Ki + R1 

to 
    f1 = k1 * u1(t) * S2

the Michaelis Menten component is simplified to k1,
R1 is not included any more.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

# initial state
x0 = [1.03, .38, .62, .52, .52]

# to make u1, u2 as functions
def u1(t):
    return 1

def u2(t):
    return 1

N = np.loadtxt("pyADAPT/examples/toy/N3.txt")
C = np.diag([1, 1, 1, 1, 0])


def reactions(t, x, k1):
    S1 = x[0]
    S2 = x[1]
    S3 = x[2]
    S4 = x[3]

    # parameters
    k2 = 1
    k3 = 0.1
    k4 = 0.5
    k5 = 1

    f1 = k1 * u1(t) * S2
    f2 = k2 * u2(t) * S3
    f3 = k3 * S1
    f4 = k4 * S1
    f5 = k5 * S4

    return np.array([f1, f2, f3, f4, f5])


def toyODE(t, x, k1):
    fx = reactions(t, x, k1)
    dxdt = N.dot(fx)
    return dxdt


# optimizing parameter k1

