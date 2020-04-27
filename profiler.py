import cProfile

from pyADAPT import optimize, DataSet
from pyADAPT.examples.toy import ToyModel

if __name__ == '__main__':
    toy = ToyModel()
    data = DataSet(raw_data_path='data/toyModel/toyData.mat',
                   data_specs_path='data/toyModel/toyData.yaml')
    # print(isinstance(toy, Model))
    optimize(toy, data, 'k1', n_core=1, n_iter=20, delta_t=0.1)

    cProfile.run("optimize(toy, data, 'k1', n_iter=20, delta_t=0.1)",
                 filename="pyADAPT.prof")
