""" Data set format specification
ADAPT as a stochastic simulation method, depends heavily on the data set format.
The following dataset is good enough for the toy model, but depends heavily on
the MATLAB data file format.

Current, xarray/pandas and pickle seems promising

The dataset provided to ADAPT procedure should always be a list of distributions.
But the storing and exchanging format should be pandas pickle.

This is not an urgent problem to solve, I should first try to use hard coded
script to provide the data set structures.

TODO: routine `read_pandas` (read_data)
"""

import pickle

import numpy as np

from pyADAPT.state import State
from pyADAPT.io import read_data_raw, read_data_specs


class DataSet(list):
    """
    dataset for ADAPT. containing phenotypes measured at different stages after
    intervention. list of States

    DataSet only stores the source data and yields new interpolant data. Model
    stores the interpolants as the trajectories.
    """
    def __init__(self,
                 raw_data_path="",
                 data_specs_path="",
                 raw_data={},
                 data_specs={},
                 name=""):
        """
        raw_data: phenotypes organized into a dictionary

        data_info: instructions of the data, such as which time variable
            should be used for which state.
        """
        self.data_specs = data_specs if data_specs else read_data_specs(
            data_specs_path)
        self.name = name if name else self.data_specs["groups"][0]
        self.raw_data = (raw_data[self.name] if raw_data else read_data_raw(
            raw_data_path, group=self.name))

        self.structure = {}
        self.ordered_names = []

        assert "structure" in self.data_specs

        for k, v in self.data_specs.items():
            self.__setattr__(k, v)

        for k, v in self.structure.items():
            time = self.raw_data[v["time"]]
            means = self.raw_data[v["means"]]
            stds = self.raw_data[v["stds"]]

            try:
                time_unit = v["time_unit"]
            except KeyError as e:
                time_unit = "seconds"
                print(
                    f"Warning: undefined {e.args[0]}, fallback to default ({time_unit})"
                )

            try:
                unit = v["unit"]
            except KeyError as e:
                unit = "mM"
                print(
                    f"Warning: undefined {e.args[0]}, fallback to default ({unit})"
                )

            s = State(
                name=k,
                time=time,
                time_unit=time_unit,
                means=means,
                stds=stds,
                unit=unit,
            )

            self.append(s)
            self.ordered_names.append(k)

    def get_time(self, when):
        if when == "end":
            time = min([s.time[-1] for s in self])
        elif when == "begin":
            time = max([s.time[0] for s in self])
        return time

    @property
    def end_time(self):
        return self.get_time("end")

    @property
    def begin_time(self):
        return self.get_time("begin")

    def get_state_names(self):
        return [s.name for s in self]

    def interpolate(self, n_ts=100, method="Hermite") -> np.ndarray:
        """In every ADAPT iteration, this function is called once to get a new
        spline for the optimizer to fit (from t0 till the end). the length of
        the list of the splines should equal the number of states in the data.

        return
        ------
        numpy.ndarray in the same order as `self.ordered_names`
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
        inter_p = np.zeros((len(self), n_ts, 2))  # a stands for array
        # inter_p [0,:, 0] = np.ones(n_ts) * self
        # v = list(self.values())
        for i in range(len(self)):
            inter_p[i, :, 0] = self[i].interp_values(n_ts=n_ts)
            inter_p[i, :, 1] = self[i].interp_stds(n_ts=n_ts)
        return inter_p

    def __getitem__(self, index):
        if type(index) is str:
            index = self.ordered_names.index(index)
        return super().__getitem__(index)


def plot_splines(D, N, n_ts=100, axes=None, seed=0):
    if axes is not None:
        # assert axes.size == len(D)
        ncols = axes.shape[1]
        fig = axes.flatten()[0].get_figure()
        plt = None
    else:
        import matplotlib.pyplot as plt

        plt.style.use("ggplot")
        ncols = get_cols(len(D))
        nrows = int(np.ceil(len(D) / ncols))
        fig, axes = plt.subplots(nrows, ncols, squeeze=False)
        fig.canvas.set_window_title(f"Interpolation of {D.name}")

    np.random.seed(seed)
    for i in range(N):
        idp = D.interpolate(n_ts=n_ts)
        for j in range(idp.shape[0]):
            ax = axes[j // ncols, j % ncols]
            state = D[j]
            t_ = state.time
            t = np.linspace(t_[0], t_[-1], n_ts)
            ax.plot(t, idp[j, :, 0], color="red", alpha=0.15)
            state.plot_samples(ax)

            if not ax.title.get_label():
                state.errorbar(ax)
    fig.tight_layout()
    if plt:
        plt.show()


def get_cols(N, ratio=1):
    return int(np.ceil(np.sqrt(N * ratio)))


if __name__ == "__main__":
    from pprint import pprint, pformat

    # D = DataSet(
    #     raw_data_path="data/toyModel/toyData.mat",
    #     data_specs_path="data/toyModel/toyData.yaml",
    # )
    D = DataSet(
        raw_data_path="data/trehalose/smallbone2011_data.mat",
        data_specs_path="data/trehalose/smallbone2011_data.yaml",
    )
    idp = D.interpolate(n_ts=10)
    # all for s1:
    pprint(idp[0, :, :])
    # select all the data from time step 3
    pprint(idp[:, 3, :])
    # all for values
    pprint(idp[:, :, 0])
    # all for stds
    pprint(idp[:, :, 1])

    n_interp = 100
    n_ts = 200
    plot_splines(D, n_interp, n_ts)
