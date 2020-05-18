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
#%%
from scipy.optimize import least_squares, leastsq
import numpy as np
import pandas as pd
from pyADAPT.examples import Smallbone2011
from van_heerden_preprocess import vhd_dataset


smallbone = Smallbone2011()

n_ts = 100
data = dict()
for s in vhd_dataset:
    data[s.name] = s.interp_means(n_ts)

#%% plot data
import matplotlib.pyplot as plt
plt.style.use(["science", "grid"])
fig, axes = plt.subplots(2, 3, figsize=(9, 6))
# with plt.style.context("science"):
for i, s in enumerate(vhd_dataset.names):
    ax:plt.Axes = axes.flatten()[i]
    ax.plot(vhd_dataset[s].get_timepoints(n_ts), data[s])
    ax.set_title(s)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('concentration (mM)')

fig.tight_layout()


# %%
