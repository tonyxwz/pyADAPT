import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numpy.random import rand, randn
from ode_true import ode_true
from fluxes_true import fluxes_true


def main():
    plt.figure()
    im = plt.imread('true_system.png')
    plt.imshow(im)
    x0 = [1.03, .38, .62, .52, .52]
    tspan = [-1e3, 0]

    k6range = np.logspace(-2, -3.5, 5)
    k60 = k6range[0]
    
    ydata = []
    v = []
    for k6 in k6range:
        sol = solve_ivp(lambda t, x: ode_true(t, x, k6),
                    tspan, x0)
        # print(sol.y[:,-1])
        ydata.append(sol.y[:,-1])
        v.append(fluxes_true(0, ydata[-1]))
    
    ydata = np.array(ydata)
    v = np.array(v)
    
    plt.figure()
    plt.plot(range(5), v, "*-")
    plt.legend(["v1", "v2", "v3", "v4", "v5"])
    plt.title('fluxes')
    
    
    stddata = .2 * np.diag(5)


if __name__ == "__main__":
    main()