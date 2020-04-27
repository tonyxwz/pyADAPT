""" manipulating 3-D data with coordinates: iter, time id
some wrappers over tedious xarray functions
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


def mean(traj: xr.DataArray):
    return traj.mean(dim="iter")


def time_points(traj):
    return traj.coords['time']


def plot_mean(traj: xr.DataArray, axes=None, **plot_options):
    if axes is not None:
        axes = np.asarray(axes).flatten()
        assert len(axes) == len(traj.coords['id'])
    else:
        pass
    m = mean(traj)

    print(m)
    for i, p in enumerate(traj.coords['id']):
        axes[i].plot(time_points(traj), m[:, i], **plot_options)
        axes[i].set_title(p.data[()])


def plot(traj: xr.DataArray, axes=None, **plot_options):
    if axes is not None:
        axes = np.asarray(axes).flatten()
        assert len(axes) == len(traj.coords['id'])
    else:
        pass

    for j, p in enumerate(traj.coords['id']):
        for i_iter in traj.coords['iter']:
            axes[j].plot(time_points(traj), traj.loc[i_iter, :, p],
                         **plot_options)
