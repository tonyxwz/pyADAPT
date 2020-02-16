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
import scipy.io as scio

__all__ = ["read_data_specs", "read_data_raw", "read_mat", "save_npy"]


def read_data_specs(path_str):
    """ yaml or json """
    with open(path_str, "r") as f:
        ext = path_str.split(".")[-1]
        if ext == "yaml":
            info = load(f, Loader=Loader)
        elif ext == "json":
            info = json.load(f)
    return info


def read_data_raw(path_str, group=""):
    file_ext = path_str.split(".")[-1]
    if file_ext == "mat":
        data_dict = read_mat(matpath=path_str)
        if group:
            data_dict = data_dict[group]
    elif file_ext == "npy":
        data_dict = np.load(path_str)
    elif file_ext == "json":
        pass
    elif file_ext == "csv":
        pass
    else:
        raise Exception("Unknown file type.")

    return data_dict


def read_mat(matpath=""):
    mat = scio.loadmat(matpath, squeeze_me=False)

    del mat["__header__"]
    del mat["__version__"]
    del mat["__globals__"]

    mat2 = dict()  # {datasetname: dataset}
    for k, v in mat.items():
        mat2[k] = dict()  # {fieldName: data}
        for f in v.dtype.fields:
            mat2[k][f] = np.array([i[0] for i in v[f][0][0]])
    return mat2


def save_npy(data):
    # np.save is simple enough; no need for this function
    pass


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    mat = read_data_raw(path_str="data/clampModelFinal/DataFinal2.mat")
    for datasetname, dataset in mat.items():
        keys = dataset.keys()
        names = [k[:-5] for k in keys if (k[-4:] == "time") and (len(dataset[k]) > 1)]
        for n in names:
            try:
                time, flux = dataset[n + "_time"], dataset[n + "_flux"]
                fig, axes = plt.subplots()
                axes.set_title(n)
                axes.plot(time, flux)
            except KeyError as ke:
                print(ke)
    plt.show()
