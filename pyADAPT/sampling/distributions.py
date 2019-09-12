import numpy as np


class Normal_Dist(object):
    """ Normal distribution """
    # one distribution from the data
    def __init__(self, mean, std):
        self.mean = mean
        self.std = std

    def sample(self):
        return np.random.normal(self.mean, self.std)