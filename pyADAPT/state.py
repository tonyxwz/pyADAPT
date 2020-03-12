# -*- coding: utf-8 -*-
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline, CubicSpline
from cached_property import cached_property

from pyADAPT.sampling import NormalDist

__all__ = ["State"]

spline_map = {"Hermite": PchipInterpolator, "Cubic": CubicHermiteSpline}


class State(list):
    """ species, x, implemented as list of NormalDist
    ```
    [
        (time[0], (mean=2, std=0.7)),
        (time[1], (mean=3, std=1.5)),
        (time[2], (mean=5, std=1.1)),
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
    ):
        super().__init__()
        self.name = name
        self.time = np.array(time)
        self.time_unit = time_unit
        self.unit = unit
        self.sampled_values = None
        assert len(time) == len(stds) == len(means)
        for m, s in zip(means, stds):
            self.append(NormalDist(m, s))

    @property
    def values_spline(self):
        """ interpolate the values """
        self.sampled_values = self.sample()
        # remove this line if the first interpolation should be random
        # self.sampled_values[0] = self[0].mean
        pp = PchipInterpolator(self.time, self.sampled_values)
        return pp

    def interp_values(self, n_ts=100, method="pchip"):
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        values_interp = self.values_spline(tnew)
        return values_interp

    #  doesn't change between iterations at all
    @cached_property
    def variances_spline(self):
        """ interpolate standard deviation
        to calculate the objective (error) function
        interpolate linearly on variance
        """
        # variances = np.asarray([d.std for d in self]) ** 2
        pp = PchipInterpolator(self.time, self.variances)
        return pp

    def interp_stds(self, n_ts=100, method="pchip"):
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        variances_interp = self.variances_spline(tnew)
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
        return self.stds**2

    def sample(self):
        """ random sample of all the data points as an array """
        return np.asarray([d.sample() for d in self])

    def __repr__(self):
        repr_list = list(zip(self.time, self))
        return f"State({self.name}): {repr_list}"

    @cached_property
    def plt(self):
        import matplotlib.pyplot as plt

        plt.style.use("ggplot")
        return plt

    def errorbar(self, axes=None):
        if axes is None:
            import matplotlib.pyplot as plt

            plt.style.use("ggplot")
            _, axes = plt.subplots()
        else:
            plt = None

        axes.set_title(self.name)
        axes.set_xlabel(self.time_unit + "(s)")
        axes.set_ylabel(self.unit)
        axes.errorbar(
            self.time[1:],
            [d.mean for d in self[1:]],
            yerr=[d.std for d in self[1:]],
            fmt=".b",
            uplims=True,
            lolims=True,
            elinewidth=1,
        )
        if plt is not None:
            axes.bar(self.time, [d.mean for d in self])
        return axes

    def plot_samples(self, axes=None):
        if axes is None:
            import matplotlib.pyplot as plt

            plt.style.use("ggplot")
            _, axes = plt.subplots()
        else:
            plt = None
        axes.plot(self.time,
                  self.sampled_values,
                  ".g",
                  alpha=0.5,
                  markersize=5)

        return axes


if __name__ == "__main__":
    from pprint import pprint
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots()

    time = [0, 3, 4, 7, 10]
    means = [6, 4, 2, 5, 3]
    stds = [1, 0.5, 1.7, 1.2, 0.9]

    time_unit = "day"
    unit = "mM/L"

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
    s.plot_samples(axes=axes)

    plt.show()
