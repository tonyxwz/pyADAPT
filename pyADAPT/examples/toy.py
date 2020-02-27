import numpy as np
import pandas as pd

from pyADAPT.basemodel import BaseModel


class ToyModel(BaseModel):
    def __init__(self):
        self.name = "Toy Model"
        self.description = "Toy Model as appeared in the 2013 paper by Natal van Riel"

        self.add_parameter(name="k1", value=1, vary=True, lb=0)
        self.add_parameter(name="k2", value=1, vary=False, lb=0)
        self.add_parameter(name="k3", value=0.1, vary=False, lb=0)
        self.add_parameter(name="k4", value=0.5, vary=False, lb=0)
        self.add_parameter(name="k5", value=1, vary=False, lb=0)
        self.init_vary = ["k1", "k2", "k3", "k4", "k5"]

        self.add_state(name="s1", value=1)
        self.add_state(name="s2", value=1)
        self.add_state(name="s3", value=1)
        self.add_state(name="s4", value=1)

        super().__init__()

    def odefunc(self, t, x, p):
        """ODE function of the toy model
        `t`: time
        `x`: current state
        `p`: parameters (`self.parameters)
        `u`: constants (`self.constants`)
        """
        v1, v2, v3, v4, v5, ds3dt = self.fluxes(t, x, p)

        dxdt = np.zeros(len(x))
        dxdt[0] = v1 - v3 - v4
        dxdt[1] = -v1 + v2
        dxdt[2] = ds3dt
        dxdt[3] = v4 - v5

        return dxdt

    def fluxes(self, t, x, p):
        u = self.inputs(t)
        u1 = u["u1"]
        u2 = u["u2"]

        s1 = x[0]
        s2 = x[1]
        s3 = x[2]
        s4 = x[3]

        k1 = p["k1"]
        k2 = p["k2"]
        k3 = p["k3"]
        k4 = p["k4"]
        k5 = p["k5"]

        v = np.zeros(6)
        v[0] = k1 * u1 * s2
        v[1] = k2 * u2 * s3
        v[2] = k3 * s1
        v[3] = k4 * s1
        v[4] = k5 * s4
        v[5] = v[0] - v[1]

        return v

    def inputs(self, t):
        # let's skip `t` for the ToyModel
        return {"u1": 1, "u2": 1}


if __name__ == "__main__":
    import os
    from pyADAPT.dataset import DataSet
    from pyADAPT.optimize import optimize
    import cProfile

    model = ToyModel()
    data = DataSet(
        raw_data_path="data/toyModel/toyData.mat",
        data_specs_path="data/toyModel/toyData.yaml",
    )
    # cProfile.run("optimize(model, data, 'k1', n_iter=1, n_tstep=10)")
    # d = data.interpolate()
    # print("begin")
    # k = model.compute_states(0.2, [0, 2], d[:, 0, 0], new_param_names=["k1"])
    # print(k)
    ptraj, straj = optimize(model, data, "k1", n_iter=10, n_tstep=100)
