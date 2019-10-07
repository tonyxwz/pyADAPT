import unittest

import numpy as np

from pyADAPT import DataSet
from pyADAPT.io import read_data_raw


class Test_Read_npy(unittest.TestCase):
    def test_unittest(self):
        self.assertNotEqual(20, 'ere')

    def test_read(self):
        npy_path = 'data/toyModel/toyData.npy'
        result = read_data_raw(npy_path)
        self.assertIsInstance(result, np.ndarray)

    def test_interpolate(self):
        n_ts = 187
        D = DataSet(raw_data_path='data/toyModel/toyData.npy',
                    data_info_path='data/toyModel/toyData.yaml')
        spl = D.interpolate(n_ts=n_ts)
        self.assertEqual(len(spl), 4)
        s1 = spl['s1']
        self.assertEqual(len(s1['values']), n_ts)
        self.assertEqual(len(s1['stds']), n_ts)
