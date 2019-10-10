import numpy as np

from pyADAPT.core import State
from pyADAPT.io import read_data_info, read_data_raw


class DataSet(dict):
    """
    dataset for ADAPT. containing phenotypes measured at different stages after
    intervention. list of States

    Job division from DataSet and Model: DataSet only stores the source data and
    yields new interpolant data. Model stores the interpolants as the
    trajectories.
    """
    def __init__(self, raw_data_path="", data_info_path="",
                 raw_data={}, data_info={}):
        """
        raw_data: phenotypes organized into a dictionary

        data_info: instructions of the data, such as which time variable 
            should be used for which state.
        """
        if data_info:
            self.info = data_info
        else:
            self.info = read_data_info(data_info_path)

        if raw_data:
            self.data = raw_data
        else:
            # group: to solve the stupid matlab issue
            group = self.info['mat_group'] if 'mat_group' in self.info else ""
            self.data = read_data_raw(raw_data_path, group=group)

        self.structure = {}
        for k,v in self.info.items():
            self.__setattr__(k, v)

        for k, v in self.structure.items():
            # !print(v)
            time = self.data[v['time']]
            means = self.data[v['means']]
            stds = self.data[v['stds']]
            s = State(name=k, time=time, means=means, stds=stds)
            self[k] = s

    def interpolate(self, n_ts=100, method='Hermite'):
        """In every ADAPT iteration, this function is called once to get a new
        spline for the optimizer to fit (from t0 till the end). the length of
        the list of the splines should equal the number of states in the data.

        return
        ------
        numpy.ndarray (In py3.7+, OrderedDict is not required.)
        ```
        [
            [
                [value, std],
                ......
                [value, std]
            ]
        ]
        ```
        """
        # TODO: add different interp methods
        # TODO: take care of states using different time points, t1, t2
        ainter_p = np.zeros((len(self), n_ts, 2))  #a stands for array
        v = list(self.values())
        for i in range(len(v)):
            ainter_p[i, :, 0] = v[i].interp_values(n_ts=n_ts)
            ainter_p[i, :, 1] = v[i].interp_stds(n_ts=n_ts)
        return ainter_p


if __name__ == "__main__":

    from pprint import pprint, pformat
    # from colorama import Fore, Back, Style, init
    # init()
    D = DataSet(raw_data_path='data/toyModel/toyData.npy',
        data_info_path='data/toyModel/toyData.yaml')
    idp = D.interpolate(n_ts=10)

    # all for s1:
    pprint(idp[0, :, :])

    # select all the data from time step 3
    pprint(idp[:,3,:])

    # all for values
    pprint(idp[:, :, 0])


