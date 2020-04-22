import unittest

import numpy as np
import pandas as pd

from pyADAPT.basemodel import BaseModel
from pyADAPT.examples.lotka import LotkaVolterra
from pyADAPT.examples.toy import ToyModel


class SimpleModel(BaseModel):
    def __init__(self):
        self.add_parameter(name="k1", value=0.1, vary=False, lb=0, ub=np.inf)
        self.add_parameter(name="k2", value=1, vary=False, lb=0, ub=np.inf)
        self.add_parameter(name="k3", value=0.45, vary=True, lb=0, ub=1)

        self.add_state(name="s1", value=1.03, observable=True)
        self.add_state(name="s2", value=0.57, observable=False)
        super().__init__()

    def ode_func(self, t, x, p):
        return super().ode_func(t, x, p)

    def fluxes_func(self, t, x, p):
        return super().fluxes_func(t, x, p)

    def inputs_func(self, t):
        return super().inputs_func(t)


class TestLotka(unittest.TestCase):
    def setUp(self):
        self.model = LotkaVolterra()

    def test_compute(self):
        y = self.model.compute_states(t_span=[0, 100],
                                      x0=[10, 10],
                                      t_eval=np.linspace(0, 100, 500))
        self.assertEqual(y.shape, (2, 500))


class TestBaseModel(unittest.TestCase):
    def setUp(self):
        self.model = SimpleModel()

    def test_components(self):
        print(self.model.parameters[self.model.parameters["vary"]])
        self.assertIsInstance(self.model.parameters, pd.DataFrame)
        self.assertIsInstance(self.model.states, pd.DataFrame)
        self.assertEqual(self.model.parameters.loc["k1"]["name"], "k1")
        self.assertEqual(self.model.states.loc["s1"]["value"], 1.03)
