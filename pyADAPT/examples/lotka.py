import numpy as np
import math
from pyADAPT.basemodel import BaseModel
from pyADAPT.dataset import DataSet
from pyADAPT.optimize import Optimizer, optimize


class LotkaVolterra(BaseModel):
    """ Lotka Volterra Model using ADAPT

    https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations
    http://www.tiem.utk.edu/~gross/bioed/bealsmodules/predator-prey.html


    Predictor
    t: time

    states:
    x: the number of prey (rabbit)  dxdt = αx - βxy
    y: the number of predator (fox) dydt = δxy -  γy

    Parameters:

    α, β, γ, δ: positive real parameters describing the interaction between
    the two populations
    """
    def __init__(self):
        self.add_state(name="prey", value=10, observable=True)
        self.add_state(name="predator", value=10, observable=True)

        self.add_parameter(name="alpha", value=1.1, vary=False, lb=0)
        self.add_parameter(name="beta", value=0.4, vary=False, lb=0)
        self.add_parameter(name="delta", value=0.1, vary=False, lb=0)
        self.add_parameter(name="gamma", value=0.4, vary=False, lb=0)
        super().__init__()

    def odefunc(self, t, x, p):
        alpha = p["value"].loc["alpha"]
        beta = p["value"].loc["beta"]
        delta = p["value"].loc["delta"]
        gamma = p["value"].loc["gamma"]

        d_prey = alpha * x[0] - beta * x[0] * x[1]
        d_predator = delta * x[0] * x[1] - gamma * x[1]

        return np.array([d_prey, d_predator])


class WrongLotkaVolterra(BaseModel):
    """ If the parameters are not accurate, simply use optimization method.
    In ADAPT, the structure of the model should be wrong not just parameter
    being inaccurate, change the odefunc a bit
    """
    def __init__(self):
        self.add_state(name="prey", value=10, observable=True)
        self.add_state(name="predator", value=10, observable=True)

        self.add_parameter(name="alpha", value=1.1, vary=False, lb=0)
        self.add_parameter(name="beta", value=0.4, vary=False, lb=0)
        self.add_parameter(name="delta", value=0.1, vary=False, lb=0)
        self.add_parameter(name="gamma", value=0.4, vary=False, lb=0)
        super().__init__()

    def odefunc(self, t, x, p):
        alpha = p.at["alpha", "value"]
        beta = p.at["beta", "value"]
        delta = p.at["delta", "value"]
        gamma = p.at["gamma", "value"]

        # d_prey = alpha * x[0] - beta * x[0] * x[1]
        # note that the modification cannot be complemented by a factor
        d_prey = alpha * x[0] - beta * x[0] * x[1] * math.log(x[1] + 1)
        d_predator = delta * x[0] * x[1] - gamma * x[1]

        return np.array([d_prey, d_predator])


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    t0 = 0
    tf = 100
    n = 500

    lotka = LotkaVolterra()
    y = lotka.compute_states(t_span=[t0, tf],
                             x0=[10, 10],
                             t_eval=np.linspace(t0, tf, n))

    noise = 0.2 * np.random.randn(2, n) * y
    data = noise + y

    plt.plot(np.linspace(t0, tf, n), data[0, :])
    plt.plot(np.linspace(t0, tf, n), data[1, :], "--")
    plt.legend(lotka.states["name"])
    plt.title("True model output")

    wrong_lotka = WrongLotkaVolterra()
    y2 = wrong_lotka.compute_states(t_span=[t0, tf],
                                    x0=[10, 10],
                                    t_eval=np.linspace(t0, tf, n))
    plt.figure()
    plt.plot(np.linspace(t0, tf, n), y2[0, :])
    plt.plot(np.linspace(t0, tf, n), y2[1, :], "--")
    plt.legend(wrong_lotka.states["name"])
    plt.title("Wrong model output")

    # TODO create fake dataset, apply ADAPT here
    plt.show()
