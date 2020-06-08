"""
some utilities for preprocessing the data
"""

from scipy.ndimage import gaussian_filter1d


def smooth(data):
    # figure out gaussian order and stdev automatically from the data
    return gaussian_filter1d
