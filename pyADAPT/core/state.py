#!/usr/bin/env python
# -*- coding: utf-8 -*-
# from cobra.core import Species
from scipy import interpolate
import numpy as np


from pyADAPT.sampling import Normal_Dist


class State(list):
    """ species, x, implemented as list of Distribtions """
    def __init__(self, name="", time=None, means=None, stds=None, observable=False):
        super().__init__()
        self.name = name  # s1
        self.time = np.array(time)  # t1 from dataset
        print(time)
        print(stds)
        print(means)
        assert len(time) == len(stds) == len(means)
        for m, s in zip(means, stds):
            self.append(Normal_Dist(m, s))

        self.observable = observable  # IDK what does this do at all...
        # self.spline = self.generate_spline()
        self.value_hist = list()
        self.std_hist = list()

    def value_spline(self, n_ts=100):
        """ interpolate the values """
        values = self.sample()
        tck = interpolate.splrep(self.time, values)
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        xnew = interpolate.splev(tnew, tck, der=0)
        self.value_hist.append(values)
        return xnew

    def std_spline(self, n_ts=100):
        """ interpolate standard deviation 
        To be used to calculate the objective (error) function
        interpolate linearly on variance
        """
        variances = np.asarray([d.std for d in self]) ** 2
        tck = interpolate.splrep(self.time, variances)
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        variances_new = interpolate.splev(tnew, tck, der=0)
        std_new = np.sqrt(variances_new)
        self.std_hist.append(std_new)
        return std_new


    def sample(self):
        """ random sample of all the data points as an array """
        return np.asarray([ d.sample() for d in self ])
