#!/usr/bin/env python
# -*- coding: utf-8 -*-
# from cobra.core import Species

import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline

from pyADAPT.sampling import NormalDist


class State(list):
    """ species, x, implemented as list of NormalDist
        [
            (time[0], (mean=2, std=0.7)),
            (time[1], (mean=3, std=1.5)),
            (time[2], (mean=5, std=1.1)),
            ...
        ]
    """
    def __init__(self, name="", time=None, means=None, stds=None, observable=False):
        super().__init__()
        self.name = name  # s1
        self.time = np.array(time)  # e.g. toymodel, t1 from dataset

        assert len(time) == len(stds) == len(means)
        for m, s in zip(means, stds):
            self.append(NormalDist(m, s))
        # TODO: seems that observable is not used here
        self.observable = observable  # IDK what does this do at all...
        # self.spline = self.generate_spline()
        self.value_hist = list()
        self.std_hist = list()

    def value_spline(self, n_ts=100):
        """ interpolate the values """
        values = self.sample()
        pp = PchipInterpolator(self.time, values)
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        xnew = pp(tnew)
        self.value_hist.append(values)
        return xnew

    # TODO: consider define std splines as a property for they don't change between iterations
    #  between iterations at all
    def std_spline(self, n_ts=100):
        """ interpolate standard deviation
        To be used to calculate the objective (error) function
        interpolate linearly on variance
        """
        variances = np.asarray([d.std for d in self]) ** 2
        pp = PchipInterpolator(self.time, variances)
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        variances_new = pp(tnew)
        std_new = np.sqrt(variances_new)
        self.std_hist.append(std_new)
        return std_new


    def sample(self):
        """ random sample of all the data points as an array """
        return np.asarray([ d.sample() for d in self ])

    def __repr__(self):
        repr_list = list()
        for i in range(len(self)):
            repr_list.append((self.time[i], self[i]))

        return f"State {self.name}: {repr_list}"


if __name__ == "__main__":
    time  =  [0,   3,   4,   7,  10]
    means =  [6,   4,   2,   5,   3]
    stds  =  [1, 0.5, 1.7, 1.2, 0.9]

    s = State(name='test_state', time=time, means=means, stds=stds)
    print(s)
