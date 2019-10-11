import unittest

import numpy as np

from pyADAPT import DataSet
from pyADAPT.io import read_data_raw


class Test_DataSet(unittest.TestCase):
    def setUp(self):
        self.npy_path = 'data/toyModel/toyData.npy'
        self.info_path = 'data/toyModel/toyData.yaml'
        self.D = DataSet(raw_data_path='data/toyModel/toyData.npy',
                    data_specs_path='data/toyModel/toyData.yaml')

    def test_read(self):
        result = read_data_raw(self.npy_path)
        self.assertIsInstance(result, np.ndarray)

    def test_interpolate(self):
        n_ts = 187
        spl = self.D.interpolate(n_ts=n_ts)
        self.assertEqual(spl.shape, (4, n_ts, 2))
