from cobra.core import Species
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
        self.history = list()

    def generate_spline(self, n_ts=100):
        values = self.sample()
        tck = interpolate.splrep(self.time, values)
        tnew = np.linspace(self.time[0], self.time[-1], n_ts)
        xnew = interpolate.splev(tnew, tck, der=0)
        self.history.append(values)
        return xnew

    def sample(self):
        # random sample of all the data points as an array
        v = np.zeros(len(self))
        for i in range(len(self)):
            v[i] = self[i].sample()
        return v
