import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from pyADAPT.examples.toy.case1 import toyODE
from pyADAPT.examples.toy.case1 import x0

# noise:
# np.random.randn(size)

def noise(size):
    ub = np.array(x0 * 0.2)
    lb = np.zeros(len(x0))
    sigmas = [np.random.rand() * x for x in (ub - lb)]
    xi = [ np.random.randn() * x for x in sigmas ]
    return xi


def generate_data(x0):
    fig, axes = plt.subplots(2, 2)
    for i in np.arange(5):
        print(x0)
        t_span = [i*1000, (i+1) * 1000]
        sol = solve_ivp(lambda t, y: toyODE(t, y, i), t_span, x0)
        # print(sol.y.shape)
        for j in range(4):
            ax = axes[j//2, j%2]
            if not ax.title.get_label():
                ax.set_title('s' + str(j+1))
            ax.plot(sol.t, sol.y[j, :])
        x0 = sol.y[:, -1]
    fig.tight_layout()

if __name__ == "__main__":
    generate_data(x0)
    plt.show()
