import xarray as xr


class Trajectory(xr.DataArray):
    def __init__(self, **kw):
        super().__init__(**kw)

    def plot(self):
        pass

    def means(self):
        pass
