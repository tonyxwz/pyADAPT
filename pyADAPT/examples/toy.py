import numpy as np
import pandas as pd

from pyADAPT.basemodel import BaseModel


class ToyModel(BaseModel):
    def __init__(self):
        self.name = "Toy Model"
        self.notes = "Toy Model as appeared in the 2013 paper by Natal van Riel"

        self.add_parameter(name="k1", value=1, vary=True, lb=0)
        self.add_parameter(name="k2", value=1, vary=False, lb=0)
        self.add_parameter(name="k3", value=0.1, vary=False, lb=0)
        self.add_parameter(name="k4", value=0.5, vary=False, lb=0)
        self.add_parameter(name="k5", value=1, vary=False, lb=0)

        super().__init__(state_order=['s1', 's2', 's3', 's4'],
                         flux_order=['v1', 'v2', 'v3', 'v4', 'v5', 'v6'],
                         input_order=['u1', 'u2'])

    def state_ode(self, t, x, p):
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
        u2 = u[self['u2']]
        # also possible if don't mind the "no-member" lint error
        u1 = u[self.u1]

        s1 = x[self['s1']]
        s2 = x[self['s2']]
        s3 = x[self['s3']]
        s4 = x[self['s4']]

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
        return [1, 1]


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

    ptraj, straj, time = optimize(model,
                                  data,
                                  "k1",
                                  n_iter=2,
                                  n_tstep=50,
                                  n_core=1)

    ptraj, straj, time = optimize(model,
                                  data,
                                  "k1",
                                  n_iter=2,
                                  n_tstep=50,
                                  n_core=4)
