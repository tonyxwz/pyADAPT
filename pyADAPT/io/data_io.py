""" It's not possible to standardize the data flow from disk to memory.
At least the process after the data is loaded should be standardize into using
numpy and python dictionaries only. 
"""
import numpy as np
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
import json


def read_data_info(path_str):
    """ yaml or json """
    with open(path_str, 'r') as f:
        ext = path_str.split('.')[-1]
        if ext == 'yaml':
            info = load(f, Loader=Loader)
        elif ext == 'json':
            info = json.load(f)
    return info


def read_data_raw(path_str, group=""):
    file_ext = path_str.split('.')[-1]
    if file_ext == 'mat':
        import matlab.engine
        eng = matlab.engine.start_matlab()
        mat_data_dict = eng.load(path_str, nargout=1)
        if len(group):
            mat_data_dict = mat_data_dict[group]
        data_dict = dict()
        for k,v in mat_data_dict.items():
            data_dict[k] = np.squeeze(np.asarray(v))
    elif file_ext == 'npy':
        data_dict = np.load(path_str)
    elif file_ext == 'json':
        pass
    elif file_ext == 'csv':
        pass
    else:
        raise Exception('Unknown file type.')

    return data_dict