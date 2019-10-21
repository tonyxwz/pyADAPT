import matplotlib.pyplot as plt
import numpy as np

from pyADAPT import ADAPT, DataSet
from pyADAPT.models import ToyModel
from logging import DEBUG

if __name__ == '__main__':
    toy = ToyModel()
    data = DataSet(raw_data_path='data/toyModel/toyData.mat',
                data_specs_path='data/toyModel/toyData.yaml')
    # print(isinstance(toy, Model))
    adapt = ADAPT(toy, data)
    adapt.set_options(n_iter=20, n_ts=100)
    adapt.run(n_core=4)
