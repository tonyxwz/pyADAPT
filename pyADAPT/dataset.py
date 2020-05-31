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
from pyADAPT.visualize import check_axes_single


class DataSet(list):
    """
    dataset for ADAPT. containing phenotypes measured at different stages after
    intervention. list of States

    Job division from DataSet and Model: DataSet only stores the source data and
    yields new interpolant data. Model stores the interpolants as the
    trajectories.
    """

    def __init__(
        self,
        raw_data_path="",
        data_specs_path="",
        raw_data={},
        data_specs={},
        name="",
        padding=True,
    ):
        """
        - raw_data: phenotypes organized into a dictionary
        - data_info: instructions of the data, such as which time variable
            should be used for which state.
        - padding: padding before t0, default no padding (trust the data)
        """
        self.data_specs = data_specs if data_specs else read_data_specs(data_specs_path)

        #  NOTE: there're a lot to improve on this messy init method, but it works
        self.name = name if name else self.data_specs["groups"][0]
        self.raw_data = (
            raw_data[self.name]
            if raw_data
            else read_data_raw(raw_data_path, group=self.name)
        )
        self.padding = padding
        self.structure = self.data_specs["structure"]

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
                padding=self.padding,  # propagate to States
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
        # FIXME assert that all the first knot of all the splines are the same
        return max([s.time[0] for s in self])

    def get_timepoints(self, n_ts):
        return np.linspace(self.begin_time, self.end_time, n_ts)

    def get_steplength(self, n_ts):
        return (self.end_time - self.begin_time) / n_ts

    def interpolate(self, n_ts=100, method="pchip") -> np.ndarray:
        """In every ADAPT iteration, this function is called once to get a new
        spline for the optimizer to fit (from t0 till the end). the length of
        the list of the splines should equal the number of states in the data.

        only supporting pchip interpolator now because it preserves monotonicity
        # TODO how to guarantee f'(t0) = 0; t0 != time[0]
        return
        ------
        numpy.ndarray in the same order as `self`
        """
        inter_p = np.zeros((len(self), n_ts, 2))
        # FIXME the interps below should use `time`
        time = self.get_timepoints(n_ts)
        for i in range(len(self)):
            # inter_p[i, :, 0] = self[i].interp_values(n_ts=n_ts)
            # inter_p[i, :, 1] = self[i].interp_stds(n_ts=n_ts)
            inter_p[i, :, 0] = self[i].value_spline(time)
            inter_p[i, :, 1] = self[i].std_spline(time)
        return inter_p

    def __getitem__(self, index):
        if type(index) is str:
            index = self.names.index(index)
        return super().__getitem__(index)

    def fancy_plot(self, n_samples=None, n_ts=None, axes=None):
        # plot splines
        # plot errorbar
        # emphasis samples
        self.plot_splines(n_samples=n_samples, n_ts=n_ts, axes=axes)
        for i, s in enumerate(self):
            ax = axes[i]
            s.errorbar(ax)

    def plot_splines(self, n_ts=100, n_samples=100, axes=None):
        time = self.get_timepoints(n_ts)
        ts = self.get_steplength(n_ts)
        for i in range(n_samples):
            splines = self.interpolate(n_ts=n_ts)
            for j, s in enumerate(self):
                ax: plt.Axes = axes[j]
                x = splines[j, :, 0]
                if not ax.title:
                    ax.set_title(s.name)
                ax.plot(time, x, color="red", alpha=0.15)
                spots = np.searchsorted(time, s.time, side="left")
                ax.plot(s.time, x[spots], ".g", alpha=0.15)


def get_cols(N, ratio=1):
    return int(np.ceil(np.sqrt(N * ratio)))


if __name__ == "__main__":
    from pprint import pprint, pformat

    D = DataSet(
        raw_data_path="../data/toyModel/toyData.mat",
        data_specs_path="../data/toyModel/toyData.yaml",
        padding=True,
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

    order = ["s2", "s4", "s1", "s3"]
    D.align(order)
