import numpy as np
from fluxes_true import fluxes_true


def ode_true(t, x, k6):
    v = fluxes_true(t, x)
    # k6 = p
    kd = 0.01

    N = np.array([[1, 0, -1, -1, 0], [-1, 1, 0, 0, 0], [1, -1, 0, 0, 0],
                  [0, 0, 0, 1, -1]])

    dxdt = np.zeros(5)
    dxdt[0:4] = N.dot(v)
    dxdt[4] = k6 * x[3] - kd * x[4]

    return dxdt

if __name__ == "__main__":
    print(ode_true(1, [1.03, .38, .62, .52, .52], 0.04))

