import unittest
import numpy as np
from pyADAPT import Model, DataSet, ADAPT
from examples import ToyModel


class Test_ADAPT(unittest.TestCase):
    def setUp(self):
        model = ToyModel()
        data = DataSet(raw_data_path='data/toyModel/toyData.npy',
            data_specs_path='data/toyModel/toyData.yaml')
        self.app = ADAPT(model, data)

    def test_set_options(self):
        self.assertRaises(Exception, self.app.set_options, hh=120)
