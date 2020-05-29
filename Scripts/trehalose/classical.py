# -*- coding: utf-8 -*-
"""
=====================================
classic approach parameter estimation
=====================================

Classical parameter estimation should estimate a confidence interval for each
parameter. In the Smallbone model, all the "Vmax" parameters. This can be done
in Monte Carlo's fashion, to run the simulation many times with truely random
initial states and parameters. Finally, construct the confidence intervals from
these results.

fit parameters:
    ['pgi_Vmax', 'hxt_Vmax', 'hxk_Vmax', 'pgm_Vmax', 'tpp_Vmax',
        'tps_Vmax', 'nth_Vmax', 'ugp_Vmax']
states for comparison: ['g6p', 'g1p', 'udg', 't6p', 'trh']

unlike the ADAPT method, classical approach don't need the random splines, but
splines interpolated directly from van Heerden's dataset.
"""
# %%
from pyADAPT.examples import Smallbone2011
from van_heerden_preprocess import vhd
import matplotlib.pyplot as plt
import numpy as np

vhd_dataset = vhd()
smallbone = Smallbone2011()

n_ts = 100
data = dict()
time_points = np.linspace(vhd_dataset.begin_time, vhd_dataset.end_time, n_ts)
for s in vhd_dataset:
    data[s.name] = s.mean_spline(time_points)

# %% plot data
plt.style.use(["science", "grid"])
fig = plt.figure(figsize=(9, 6))
ax1 = plt.subplot2grid((2, 6), loc=(0, 1), colspan=2)
ax2 = plt.subplot2grid((2, 6), (0, 3), colspan=2)
ax3 = plt.subplot2grid((2, 6), (1, 0), colspan=2)
ax4 = plt.subplot2grid((2, 6), (1, 2), colspan=2)
ax5 = plt.subplot2grid((2, 6), (1, 4), colspan=2)
axes = [ax1, ax2, ax3, ax4, ax5]

for i, s in enumerate(vhd_dataset.names):
    ax: plt.Axes = axes[i]
    ax.plot(time_points, data[s])
    ax.set_title(s)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("concentration (mM)")

fig.tight_layout()
fig.savefig("classical.png", dpi=200)

# %%
from pyADAPT.mc import simulate, steady_states
from scipy.optimize import least_squares
from pyADAPT.optimize import optimize, Optimizer


fit_params = list()
for id in smallbone.parameters.index:
    if id[-5:] == "_Vmax":
        fit_params.append(id)
print(fit_params)
opt = Optimizer(smallbone, vhd_dataset, parameter_names=fit_params)

parameters = opt.parameters

data2 = list()
for s in smallbone.state_order:
    if s in data:
        data2.append(data[s])
data2 = np.array(data2)
# %%
# Solve the state_ode in smallbone model on the entire time range
x0 = np.concatenate([[smallbone.initial_states[0]], data2[:, 0]])
smallbone.parameters.loc[
    [
        "pgi_Vmax",
        "hxt_Vmax",
        "hxk_Vmax",
        "pgm_Vmax",
        "tpp_Vmax",
        "tps_Vmax",
        "nth_Vmax",
        "ugp_Vmax",
    ],
    "value",
] = [13.4667, 3.67, 4.75, 100, 81.45, 1000, 100, 36.8200]
y = smallbone.compute_states(time_points=time_points, x0=x0,)

# %% plot
# fig: plt.Figure = plt.figure(figsize=(9, 6))
# gs = fig.add_gridspec(2, 6)
# fig.add_subplot(gs[0, 0:2])
# fig.add_subplot(gs[0, 2:4])
# fig.add_subplot(gs[0, 4:6])
# fig.add_subplot(gs[1, 1:3])
# fig.add_subplot(gs[1, 3:5])
# axes = fig.get_axes()

# for a
fig2, axes2 = plt.subplots(2, 3, figsize=(9, 6))
for a, y_, n in zip(axes2.flatten(), y, smallbone.state_order):
    a.plot(time_points, y_)
    a.set_title(n)
fig2.savefig("steady_states-david.png", dpi=200)

#%% calculate optimize the compute_states procedure to data2

observable = [x in data for x in smallbone.state_order]


def objective(params):
    y = smallbone.compute_states(
        new_param_names=fit_params, new_params=params, x0=x0, time_points=time_points
    )
    y[observable, :] - data2
