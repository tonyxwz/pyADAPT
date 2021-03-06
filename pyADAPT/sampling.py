# -*- coding: utf-8 -*-
import numpy as np

__all__ = ['NormalDist']

class NormalDist(object):
    """ Normal distribution """
    # one distribution from the data
    def __init__(self, mean, std):
        self.mean = mean
        self.std = std

    def sample(self):
        return np.random.normal(self.mean, self.std)

    def __repr__(self):
        return f"NormalDist(mean: {self.mean}, std: {self.std})"



if __name__ == "__main__":
    n = NormalDist(5, 1)
    print(n)