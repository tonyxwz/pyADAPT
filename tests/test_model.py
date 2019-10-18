#!C:\InstallDir\miniconda3\envs\scipy36\python.exe
"""constant parameter simulation (step 1)
fit the model to the data using constant parameters
"""
import numpy as np

from pyADAPT import Model, DataSet, ADAPT


class ToyModel(Model):
    def __init__(self):
        self.name = "Toy Model"
        self.description = "Toy Model as appeared in the 2013 paper by Natal van Riel"
        
        self.add_predictor(name='t', value=[0, 10])

        self.add_constant('u1', 1)
        self.add_constant('u2', 1)

        self.add_parameter(name='k1', value=1, vary=True)
        self.add_parameter(name='k2', value=1, vary=False)
        self.add_parameter(name='k3', value=0.1, vary=False)
        self.add_parameter(name='k4', value=0.5, vary=False)
        self.add_parameter(name='k5', value=1, vary=False)

        self.add_state(name='s1', init=1, expr="v1-v3-v4")
        self.add_state(name='s2', init=1, expr="-v1+v2")
        self.add_state(name='s3', init=1, expr="ds3dt", observable=False)
        self.add_state(name='s4', init=1, expr="v4-v5")

        self.end_extend()
        super().__init__()
    
    @staticmethod
    def odefunc(t, x, p, u):
        """ODE function of the toy model
        `t`: time
        `x`: current state
        `p`: parameters (`self.parameters)
        `u`: constants (`self.constants`)
        """
        v1, v2, v3, v4, v5, ds3dt = ToyModel.reactions(t, x, p, u)


        dxdt = np.zeros(len(x))
        dxdt[0] = v1 - v3 - v4
        dxdt[1] = -v1 + v2
        dxdt[2] = ds3dt
        dxdt[3] = v4 - v5

        return dxdt

    @staticmethod
    def reactions(t, x, p, u):
        u1 = u['u1']
        u2 = u['u2']

        s1 = x['s1']
        s2 = x['s2']
        s3 = x['s3']
        s4 = x['s4']

        k1 = p['k1']
        k2 = p['k2']
        k3 = p['k3']
        k4 = p['k4']
        k5 = p['k5']

        v = np.zeros(6)
        v[0] = k1 * u1 * s2
        v[1] = k2 * u2 * s3
        v[2] = k3 * s1
        v[3] = k4 * s1
        v[4] = k5 * s4
        v[5] = v[0] - v[1]

        return v

    def inputs(self, t):
        # let's skip this for Toy Model
        return None