""" xarray dataset """
import xarray as xr
import pandas as pd
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline


class State(pd.DataFrame):
    def __init__(self, time, *data):
        super().__init__(data)

    def interp(self, n):
        pass
