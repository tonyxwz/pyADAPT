# -*- coding: utf-8 -*-
""" Data set format specification
ADAPT as a stochastic simulation method, depends heavily on the data set format.
The following dataset is good enough for the toy model, but depends heavily on
the MATLAB data file format.

Current, xarray/pandas and pickle seems promising
"""

import pickle

import numpy as np
import matplotlib.pyplot as plt

from pyADAPT.spline import State, Flux
from pyADAPT.io import read_data_raw, read_data_specs


class DataSet(list):
    """
    dataset for ADAPT. containing phenotypes measured at different stages after
    intervention. list of States

    Job division from DataSet and Model: DataSet only stores the source data and
    yields new interpolant data. Model stores the interpolants as the
    trajectories.
    """

    def __init__(
        self, raw_data_path="", data_specs_path="", raw_data={}, data_specs={}, name=""
    ):
        """
        raw_data: phenotypes organized into a dictionary

        data_info: instructions of the data, such as which time variable
            should be used for which state.
        """
        self.data_specs = data_specs if data_specs else read_data_specs(data_specs_path)

        #  NOTE: there're a lot to improve on this messy init method, but it works
        self.name = name if name else self.data_specs["groups"][0]
        self.raw_data = (
            raw_data[self.name]
            if raw_data
            else read_data_raw(raw_data_path, group=self.name)
        )

        self.structure = {}

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
                unit = "mM/L"
                print(f"Warning: undefined {e.args[0]}, fallback to default ({unit})")
            # FIXME use Flux for flux splines
            s = State(
                name=k,
                time=time,
                time_unit=time_unit,
                means=means,
                stds=stds,
                unit=unit,
            )

            self.append(s)

    def align(self, order):
        cur = 0
        for n in order:
            if n in self.names:
                i = self.names.index(n)
                self.insert(cur, self.pop(i))
                cur += 1

    @property
    def names(self):
        return [s.name for s in self]

    @property
    def flux_names(self):
        return [s.name for s in self if type(s) is Flux]

    @property
    def state_names(self):
        return [s.name for s in self if type(s) is State]

    @property
    def end_time(self):
        return min([s.time[-1] for s in self])

    @property
    def begin_time(self):
        return max([s.time[0] for s in self])

    def get_timepoints(self, n_ts):
        return np.linspace(self.begin_time, self.end_time, n_ts)

    def interpolate(self, n_ts=100, method="Hermite") -> np.ndarray:
        """In every ADAPT iteration, this function is called once to get a new
        spline for the optimizer to fit (from t0 till the end). the length of
        the list of the splines should equal the number of states in the data.

        return
        ------
        numpy.ndarray in the same order as `self`
        ```
        [
            [value, std],,z
            ......,
            [value, std]
        ]
        ```
        """
        inter_p = np.zeros((len(self), n_ts, 2))
        # FIXME the interps below should use `time`
        time = self.get_timepoints(n_ts)
        for i in range(len(self)):
            # inter_p[i, :, 0] = self[i].interp_values(n_ts=n_ts)
            # inter_p[i, :, 1] = self[i].interp_stds(n_ts=n_ts)
            inter_p[i, :, 0] = self[i].interp_values(n_ts=n_ts)
            inter_p[i, :, 1] = self[i].interp_stds(n_ts=n_ts)
        return inter_p

    def __getitem__(self, index):
        if type(index) is str:
            index = self.names.index(index)
        return super().__getitem__(index)


def plot_splines(D, N, n_ts=100, axes=None, seed=0, figsize=(10, 10)):
    # TODO remove
    if axes is not None:
        # assert axes.size == len(D)
        fig = np.array(axes).flatten()[0].get_figure()
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
            ax = axes[j]
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
    import matplotlib.pyplot as plt

    D = DataSet(
        raw_data_path="data/toyModel/toyData.mat",
        data_specs_path="data/toyModel/toyData.yaml",
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
    fig0, axes0 = plt.subplots(2, 2, figsize=(8, 8))
    plot_splines(D, n_interp, n_ts, axes=axes0)

    order = ["s2", "s4", "s1", "s3"]
    D.align(order)
    idp = D.interpolate(n_ts=10)
    fig1, axes1 = plt.subplots(2, 2, figsize=(8, 8))
    plot_splines(D, n_interp, n_ts, axes=axes1)
    plt.show()
