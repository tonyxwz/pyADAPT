from __future__ import annotations
import xarray as xr


class Trajectory(dict):
    """Never mind..."""

    def __init__(self, data=[], coords=[], name=""):
        for i_iter in coords["iter"]:
            data = data[i_iter]

    def interp(self, t):
        pass

    def save(self):
        pass

    @staticmethod
    def load(path) -> Trajectory:
        pass
