import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from case1 import toyODE, x0
from scipy.io import loadmat, savemat


# noise:
# np.random.randn(size)


def add_noise(y_ss_i):
    # print(y_ss_i)
    ub = np.array(x0) * 0.2
    lb = np.zeros(len(x0))
    sigmas = np.array([np.random.rand() * x for x in (ub - lb)])
    xi_i = np.array([np.random.randn() * x for x in sigmas])
    return xi_i + y_ss_i


def generate_data(x_init):
    fig, axes = plt.subplots(2, 2)
    for i in np.arange(5):
        # x_init = add_noise(x_init)
        print(x_init)
        t_span = [i * 1000, (i + 1) * 1000]
        sol = solve_ivp(lambda t, y: toyODE(t, y, i), t_span, x_init)
        # print(sol.y.shape)
        for j in range(4):
            ax = axes[j // 2, j % 2]
            if not ax.title.get_label():
                ax.set_title("s" + str(j + 1))
            ax.plot(sol.t, sol.y[j, :])
        x_init = sol.y[:, -1]
    fig.tight_layout()


if __name__ == "__main__":
    from pprint import pprint
    from os.path import join, dirname

    generate_data(np.array(x0))

    toy_data = loadmat(join(dirname(__file__), "data_toy.mat"))
    # pprint(toy_data)
    concdata = toy_data["concdata"]
    fluxes = toy_data["fluxes"]
    stddata = toy_data["stddata"]
    pprint(concdata)
    pprint(stddata)
    print(concdata.shape, fluxes.shape, stddata.shape)

    fig2, axes2 = plt.subplots(2, 2)
    axes2[0, 0].bar(np.arange(5), concdata[:, 0], yerr=stddata[:5, 0])
    axes2[0, 0].set_ylim(0, 2)
    axes2[0, 1].bar(np.arange(5), concdata[:, 1], yerr=stddata[:5, 1])
    axes2[0, 1].set_ylim(0, 1)
    axes2[1, 0].bar(np.arange(5), concdata[:, 2], yerr=stddata[:5, 2])
    axes2[1, 0].set_ylim(0, 1.5)
    axes2[1, 1].bar(np.arange(5), concdata[:, 3], yerr=stddata[:5, 3])
    axes2[1, 1].set_ylim(0, 1)

    plt.show()
