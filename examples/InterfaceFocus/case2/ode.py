# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 15:10:20 2020

@author: tonyx
"""


import numpy as np
from fluxes import fluxes

def ode(t, x, k1):
    v = fluxes(t, x, k1)
    # k6 = p
    # kd = 0.01 

    N = np.array([[1, 0, -1, -1, 0], 
                  [-1, 1, 0, 0, 0], 
                  [1, -1, 0, 0, 0],
                  [0, 0, 0, 1, -1]])
    

    dxdt = N.dot(v)
    return dxdt

if __name__ == "__main__":
    print(ode(1, [1.03, .38, .62, .52, .52], 0.04))
