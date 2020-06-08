# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from cached_property import cached_property
from scipy.interpolate import CubicHermiteSpline, CubicSpline, PchipInterpolator

from pyADAPT.sampling import NormalDist
from pyADAPT.visualize import check_axes_single

__all__ = ["DataLine", "State", "Flux"]


class DataLine(list):
    """ species, x, implemented as list of NormalDist
    ```
    [
        (time[0], (mean0, std0)),
        (time[1], (mean1, std1)),
        (time[2], (mean2, std2)),
        ...
    ]
    ```
    """

    def __init__(
        self,
        name="",
        time=None,
        time_unit="second",
        means=None,
        stds=None,
        unit=None,
        observable=False,
        padding=0,
    ):
        super().__init__()
        self.name = name
        self.padding = padding
        self.time = np.array(time)
        self.time_unit = time_unit
        self.unit = unit
        self.last_sampled = None
        self.type = "DataLine"
        assert len(time) == len(stds) == len(means)
        for m, s in zip(means, stds):
            self.append(NormalDist(m, s))

    def get_timepoints(self, n_ts):
        return np.linspace(self.time[0], self.time[-1], n_ts)

    @cached_property
    def mean_spline(self):
        # This is useful for classical parameter estimation
        pp = PchipInterpolator(self.time, self.means)
        return pp

    def interp_means(self, n_ts=100):
        # classical parameter estimation
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        return self.mean_spline(tnew)

    def padded_time(self):
        if self.padding == 0:
            return self.time
        else:
            return np.concatenate(([self.time[0] - self.padding], self.time))

    @property
    def value_spline(self):
        """interpolate the values
        add a spline knot is padding is enabled such that the derivative of the
        spline at t = 0 is always 0 (pchip preserves monotonicity)
        the padding part is not evaluated
        """
        sample = self.sample()
        time = self.time
        if self.padding is True:
            time = np.concatenate(([self.time[0] - 1], time))
            sample = np.concatenate(([sample[0]], sample))
        pp = PchipInterpolator(time, sample)
        return pp

    def interp_values(self, n_ts=100, method="pchip"):
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        values_interp = self.value_spline(tnew)
        return values_interp

    #  doesn't change between iterations at all
    @cached_property
    def variance_spline(self):
        """ interpolate standard deviation
        to calculate the objective (error) function
        interpolate linearly on variance
        """
        pp = PchipInterpolator(self.time, self.variances)
        return pp

    def std_spline(self, time):
        variances_interp = self.variance_spline(time)
        stds_interp = np.sqrt(variances_interp)
        return stds_interp

    def interp_stds(self, n_ts=100, method="pchip"):
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        variances_interp = self.variance_spline(tnew)
        stds_interp = np.sqrt(variances_interp)
        return stds_interp

    @cached_property
    def means(self):
        return np.asarray([d.mean for d in self])

    @cached_property
    def stds(self):
        return np.asarray([d.std for d in self])

    @cached_property
    def variances(self):
        return self.stds ** 2

    def sample(self):
        """ random sample of all the data points as an array """
        return np.asarray([d.sample() for d in self])

    def smoothing(self, kernel):
        pass

    def errorbar(self, axes=None):
        # draw means in different color
        axes.set_title(self.name)
        axes.set_xlabel(self.time_unit)
        axes.set_ylabel(self.unit)
        axes.errorbar(
            self.time,
            [d.mean for d in self],
            yerr=[d.std for d in self],
            fmt=".b",
            uplims=True,
            lolims=True,
            elinewidth=0.3,
        )

    def plot_sample(self, axes=None):
        axes.plot(self.time, self.sample(), ".g", alpha=0.5, markersize=5)

    def __repr__(self):
        repr_list = list(zip(self.time, self))
        return f"{self.type}({self.name}): {repr_list}"


class State(DataLine):
    # ? rename to Species
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.type = "State"


class Flux(DataLine):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.type = "Flux"


if __name__ == "__main__":
    from pprint import pprint
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots()

    time = [0, 3, 4, 7, 10]
    means = [6, 4, 2, 5, 3]
    stds = [1, 0.5, 1.7, 1.2, 0.9]

    time_unit = "day"
    unit = "mM"

    s = State(
        name="Hepatic TG",
        time=time,
        means=means,
        stds=stds,
        time_unit=time_unit,
        unit=unit,
    )
    pprint(s)
    s.errorbar(axes=axes)
    _ = s.interp_values()
    s.plot_sample(axes=axes)

    plt.show()
