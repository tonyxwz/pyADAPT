import numpy as np

def fluxes_true(t, x):
    # S1, S2, S3, S4, R1
    u1 = 1
    u2 = 1
    if 5 < t < 6:
        u1 = 2

    k1 = 1
    Km = 0.1  # Km is Ki in the paper
    k2 = 1
    k3 = 0.1
    k4 = 0.5
    k5 = 1

    v1 = k1 * u1 * x[1] / (Km + x[4])
    v2 = k2 * u2 * x[2]
    v3 = k3 * x[0]
    v4 = k4 * x[0]
    v5 = k5 * x[3]

    return np.array([v1, v2, v3, v4, v5])