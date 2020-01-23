"""constant parameter simulation (step 1)
fit the model to the data using constant parameters
"""
import unittest
import numpy as np

from pyADAPT import Model, DataSet, ADAPT
from examples import ToyModel


class Test_Model(unittest.TestCase):
    def setUP(self):
        self.toy_model = ToyModel()

