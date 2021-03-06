# -*- coding: utf-8 -*-
"""
manipulating 3-D data with coordinates: iter, time id
some wrappers over tedious xarray functions
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from functools import wraps


def latexify(s):
    s = s.split("_")
    if len(s) == 1:
        return s[0]
    # it should be Vmax_{g1p} not G1p_{Vmax}
    return "$" + s[1] + "_{" + s[0] + "}" + "$"


def check_axes_traj(func):
    # `@` to create a axes is not given
    def wrapper(traj: xr.DataArray, axes=None, **kw):
        if axes is None:
            fig, axes = plt.subplots(len(traj.coords["id"]))
        else:
            fig = None
        axes = np.asarray(axes).flatten()
        assert len(axes) == len(traj.coords["id"])
        for i, x in enumerate(traj.coords["id"]):
            axes[i].set_title(latexify(x.data[()]))
        func(traj, axes=axes, **kw)
        return axes

    return wrapper


def check_axes_single(func):
    @wraps(func)
    def wrapper(*args, axes=None, **kw):
        if axes is None:
            _fig, axes = plt.subplot()
        func(*args, axes=axes, **kw)
        return axes

    return wrapper


def mean(traj: xr.DataArray):
    return traj.mean(dim="iter")


def time_points(traj):
    return traj.coords["time"].data


@check_axes_traj
def plot_mean(traj: xr.DataArray, axes=None, **plot_options):
    axes = np.asarray(axes).flatten()
    assert len(axes) == len(traj.coords["id"])
    m = mean(traj)
    for i, p in enumerate(traj.coords["id"]):
        axes[i].plot(time_points(traj), m[:, i], **plot_options)


def to_matrix(tr: xr.DataArray, size=(50, 50)):
    # convert trajectories of one variable into a matrix, 100x100 tile
    heatmap = np.zeros(size)
    # round each point in the trajectory to closest tile
    # v_range shouldn't be simply max and min
    # TODO ptraj is very light because one point is really intense
    v_range = tr.min().item(), tr.max().item()
    t_range = tr.coords["time"].data[0], tr.coords["time"].data[-1]
    time = tr.coords["time"].data

    for i_iter in tr.iter:
        line = tr.sel(iter=i_iter.item()).data
        for t, v in zip(time, line):
            x = int(
                round(((t - t_range[0]) / (t_range[1] - t_range[0])) * (size[0] - 1))
            )
            y = int(
                round(((v - v_range[0]) / (v_range[1] - v_range[0])) * (size[1] - 1))
            )
            heatmap[x, y] += 1
    return heatmap


@check_axes_traj
def plot_heatmap(
    traj: xr.DataArray,
    size=(50, 50),
    axes=None,
    cmap=plt.get_cmap("hot_r"),
    **plot_options
):
    for i, name in enumerate(traj.id):
        tr = traj.sel(id=name)
        # TODO axis ticks
        axes[i].imshow(to_matrix(tr, size=size).T, origin="lower", cmap=cmap)


@check_axes_traj
def plot_3d(traj: xr.DataArray, axes=None, **plot_options):
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111, projection="3d")
    # TODO plot frequency in 3d axes
    axes.plot_surface()


@check_axes_traj
def plot(traj: xr.DataArray, axes=None, **plot_options):
    for j, p in enumerate(traj.coords["id"]):
        for i_iter in traj.coords["iter"]:
            axes[j].plot(time_points(traj), traj.loc[i_iter, :, p], **plot_options)
