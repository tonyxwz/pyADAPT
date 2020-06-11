# -*- coding: utf-8 -*-
"""
actually run the Smallbone model on van Heerden's dataset
"""
import os
import platform
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
from mat4py import loadmat

import pyADAPT.visualize as vis
from pyADAPT.examples import Smallbone2011
from pyADAPT.io import load_traj, save_traj
from pyADAPT.optimize import ITER, Optimizer, optimize
from van_heerden_preprocess import vhd


def main():
    time_stamp = time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())
    if len(sys.argv) > 1:
        prefix = sys.argv[1] + "_"
    else:
        prefix = ""
    smallbone = Smallbone2011()

    fit_params = list()
    id: str
    for id in smallbone.parameters.index:
        params = id.split("_")
        if params[0] in ["pgi", "hxt", "hxk", "pgm", "tpp", "tps", "nth", "ugp"]:
            if params[1] not in ["shock", "activity"]:  #!
                # if id[-5:] == "_Vmax":
                fit_params.append(id)
    print(fit_params)

    # van heerden's dataset *without* padding
    vhd_dataset = vhd(padding=False)
    for para in fit_params:
        p, s, f, _t = optimize(
            smallbone,
            vhd_dataset,
            para,
            #                  alternative initial paramters from david's email (not including udg)
            #                            (Glt)                      (TPS2)   (TPS1)          (!missing)
            #                  pgi_Vmax hxt_Vmax hxk_Vmax pgm_Vmax tpp_Vmax tps_Vmax nth_Vmax ugp_Vmax
            # initial_parameters=[13.4667, 3.67, 4.75, 100, 81.45, 1000, 100, 36.8200],
            n_iter=4,  # TODO increase to a crazy number
            delta_t=2,
            #               g1p g6p trh t6p udg
            # weights=np.array([1, 1, 0.5, 1, 0.5]),
            timeout=1000,  # this should be enough for the bootstrapping
            ss_time=5000,
            max_retry=100,
            lambda_r=10,
            odesolver="LSODA",
            overture_variation="default",
        )
        p_fig, p_axes = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
        vis.plot(p, axes=p_axes, color="blue", alpha=0.2)
        vis.plot_mean(p, axes=p_axes, color="red")
        p_fig.tight_layout()
        p_fig.savefig(f"figures/p_{para}.png")

        s_fig, s_axes = plt.subplots(2, 3, figsize=(9, 5), dpi=200)
        vis.plot(s, axes=s_axes, color="blue", alpha=0.2)
        vis.plot_mean(s, axes=s_axes, color="red")
        for i, s1 in enumerate(vhd_dataset):
            ca = s_axes.flatten()[i + 1]
            s1.errorbar(ca)
        s_fig.tight_layout()
        s_fig.savefig(f"figures/states_{para}.png")

        f_fig, f_axes = plt.subplots(2, 4, figsize=(9, 5), dpi=200)
        vis.plot(f, axes=f_axes, color="blue", alpha=0.2)
        vis.plot_mean(f, axes=f_axes, color="red")
        f_fig.tight_layout()
        f_fig.savefig(f"figures/fluxes_{para}.png")

        save_traj(p, f"p_{para}_{time_stamp}.nc")
        save_traj(s, f"s_{para}_{time_stamp}.nc")
        save_traj(f, f"f_{para}_{time_stamp}.nc")


if __name__ == "__main__":
    plt.style.use(["science", "grid", "no-latex"])
    main()
