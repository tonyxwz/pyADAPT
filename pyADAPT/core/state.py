#!/usr/bin/env python
# -*- coding: utf-8 -*-
# from cobra.core import Species

import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline
from cached_property import cached_property

from pyADAPT.sampling import NormalDist

__all__ = ['State']

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
    def __init__(self, name="", time=None, means=None, stds=None, observable=False):
        super().__init__()
        self.name = name  # s1
        self.time = np.array(time)  # e.g. toymodel, t1 from dataset
        self.sampled_values = None
        assert len(time) == len(stds) == len(means)
        for m, s in zip(means, stds):
            self.append(NormalDist(m, s))

    @property  # cannot be cached because sampling is different each time
    def values_spline(self):
        """ interpolate the values """
        self.sampled_values = self.sample()
        pp = PchipInterpolator(self.time, self.sampled_values)
        return pp
    
    def interp_values(self, n_ts=100, method='pchip'):
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        values_interp = self.values_spline(tnew)
        return values_interp

    #  between iterations at all
    @cached_property
    def stds_spline(self):
        """ interpolate standard deviation
        To be used to calculate the objective (error) function
        interpolate linearly on variance
        """
        # variances = np.asarray([d.std for d in self]) ** 2
        pp = PchipInterpolator(self.time, self.variances)
        return pp

    def interp_stds(self, n_ts=100, method='pchip'):
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        variances_interp = self.stds_spline(tnew)
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
        return np.asarray([ d.sample() for d in self ])

    def __repr__(self):
        repr_list = list(zip(self.time, self))
        return f"State {self.name}: {repr_list}"


if __name__ == "__main__":
    from pprint import pprint, pformat
    time  =  [0,   3,   4,   7,  10]
    means =  [6,   4,   2,   5,   3]
    stds  =  [1, 0.5, 1.7, 1.2, 0.9]

    s = State(name='test_s1', time=time, means=means, stds=stds)
    pprint(s)
    # print(s.means, s.stds, s.variances, sep='\n')
