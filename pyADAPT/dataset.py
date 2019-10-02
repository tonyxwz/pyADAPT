import numpy as np

from .core import State
from .io import read_data_info, read_data_raw


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

    def generate_splines(self, n_ts=100):
        """In every ADAPT iteration, this function is called once to get a new
        spline for the optimizer to fit (from t0 till the end). the length of
        the list of the spline should be equal to the number of states in the
        data.
        """
        # TODO: find out influence if using b-spline rather than cubic-smooth
        spl = {}
        for k, v in self.items():
            # e.g.: k: 's1', v: s1: State
            d = dict()
            d['value'] = v.value_spline(n_ts=n_ts)
            d['std'] = v.std_spline(n_ts=n_ts)
            spl[k] = d
        return spl
