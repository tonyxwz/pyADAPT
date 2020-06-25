# -*- coding: utf-8 -*-
"""
actually run the Smallbone model on van Heerden's dataset
"""
import os
import platform
import sys
import time

import numpy as np
from mat4py import loadmat

import pyADAPT.visualize as vis

from pyADAPT.examples import Smallbone2011
from pyADAPT.optimize import ITER, Optimizer, optimize
from pyADAPT.io import load_traj, save_traj
from van_heerden_preprocess import vhd


def main():
    time_stamp = time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())
    if len(sys.argv) > 1:
        prefix = sys.argv[1] + "_"
    else:
        prefix = ""
    smallbone = Smallbone2011()
    # smallbone.parameters.loc[list(new_params.keys()), "init"] = list(
    #     new_params.values()
    # )
    # print(smallbone.parameters)
    fit_params = list()
    for id in smallbone.parameters.index:
        if id[-5:] == "_Vmax":
            fit_params.append(id)
    print(fit_params)

    # van heerden's dataset *without* padding
    vhd_dataset = vhd(padding=False)
    # for para in fit_params:
    p, s, f, _t = optimize(
        smallbone,
        vhd_dataset,
        *fit_params,
        #                  alternative initial paramters from david's email (not including udg)
        #                            (Glt)                      (TPS2)   (TPS1)          (!missing)
        #                  pgi_Vmax hxt_Vmax hxk_Vmax pgm_Vmax tpp_Vmax tps_Vmax nth_Vmax ugp_Vmax
        # initial_parameters=[13.4667, 3.67, 4.75, 100, 81.45, 1000, 100, 36.8200],
        n_iter=256,  # TODO increase to a crazy number
        delta_t=2,
        #               g1p g6p trh t6p udg
        weights=np.array([1, 1, 0.5, 1, 0.5]),
        timeout=1000,  # this should be enough for the bootstrapping
        ss_time=5000,
        max_retry=100,
        lambda_r=10,
        odesolver="Radau",
        overture_variation="lazy",
    )
    save_traj(p, f"p_{prefix}{time_stamp}.nc")
    save_traj(s, f"s_{prefix}{time_stamp}.nc")
    save_traj(f, f"f_{prefix}{time_stamp}.nc")


if __name__ == "__main__":
    main()
