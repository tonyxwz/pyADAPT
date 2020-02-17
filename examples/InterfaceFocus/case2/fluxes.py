# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 15:11:09 2020

@author: tonyx
"""


import numpy as np

def fluxes(t, x, k1):
    # S1, S2, S3, S4, R1
    
    u1 = 1
    u2 = 1
    if 5 < t < 6:
        u1 = 2

    # k1 = 1
    # Km = 0.1  # Km is Ki in the paper
    k2 = 1
    k3 = 0.1
    k4 = 0.5
    k5 = 1
    
    # k1: a lumped parameter used to capture and approximate the combined
    # processes of f1 and f6 in the true system. This parameter was
    # estimated for each phenotype independently.
    #
    # In contrust to case 1, here two sources of uncertainty were present 
    # in the inference problem. There was the noise in the data, but the
    # model was also a simplified description of reality (undermodelling).
    # **Monte Carlo** approach was used to estimate k1.

    v1 = k1 * u1 * x[1]
    v2 = k2 * u2 * x[2]
    v3 = k3 * x[0]
    v4 = k4 * x[0]
    v5 = k5 * x[3]

    return np.array([v1, v2, v3, v4, v5])
