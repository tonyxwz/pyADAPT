"""
simulate the real system

Step 1:
    generate the experimental data

"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import os

dirname = os.path.dirname(__file__)

# initial state
x0 = [1.03, .38, .62, .52, .52]

# to make u1, u2 as functions
def u1(t):
    return 1

def u2(t):
    return 1

N = np.loadtxt(os.path.join(dirname, "N.txt"))
C = np.diag([1, 1, 1, 1, 0])

def reactions(t, x, k6):

    S1 = x[0]
    S2 = x[1]
    S3 = x[2]
    S4 = x[3]
    R1 = x[4]

    # parameters
    Vmax = 1
    Ki = 0.1
    k2 = 1
    k3 = 0.1
    k4 = 0.5
    k5 = 1

    kd = 0.01


    f1 = (Vmax / (Ki + R1)) * u1(t) * S2
    f2 = k2 * u2(t) * S3
    f3 = k3 * S1
    f4 = k4 * S1
    f5 = k5 * S4
    f6 = k6 * S4 - kd * R1

    return [f1, f2, f3, f4, f5, f6]


def toyODE(t, x, i):
    k6_array = np.array([10, 4.2, 1.8, 0.7, 0.3]) * 1e-3
    k6 = k6_array[i]
    fx = reactions(t, x, k6)
    dxdt = N.dot(fx)
    return dxdt
