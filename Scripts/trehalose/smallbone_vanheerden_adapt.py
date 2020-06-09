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
    matpath = os.path.join(
        os.path.dirname(__file__),
        os.pardir,
        os.pardir,
        r"data/trehalose/parametersDavid.mat",
    )
    time_stamp = time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())
    new_params = loadmat(matpath)["p"]
    if len(sys.argv) > 1:
        prefix = sys.argv[1] + "_"
    else:
        prefix = ""
    # FIXME I think I will not use David's parameters
    # because only smallbone's parameters can reach to a steady states
    # plus, the new initialization method provids a workaround.
    rename_map = {
        "GLT_KeqGLT": None,
        "GLT_KmGLC": "hxt_Kglc",
        "GLT_VmGLT": "hxt_Vmax",
        "HXK1_Kadp": "hxk_Kadp",
        "HXK1_Katp": "hxk_Katp",
        "HXK1_Keq": "hxk_Keq",
        "HXK1_Kg6p": "hxk_Kg6p",
        "HXK1_Kglc": "hxk_Kglc",
        "HXK1_Kt6p": "hxk_Kit6p",
        "PGI1_Keq": "pgi_Keq",
        "PGI1_Kf6p": "pgi_Kf6p",
        "PGI1_Kg6p": "pgi_Kg6p",
        "PGM1_Keq": "pgm_Keq",
        "PGM1_Kg1p": "pgm_Kg1p",
        "PGM1_Kg6p": "pgm_Kg6p",
        "TPS1_Kg6p": "tps_Kg6p",
        "TPS1_Kudp_glc": "tps_Kudg",
        "TPS1_Kpi": None,
        "TPS1_KmF6P": None,
        "TPS2_Kt6p": "tpp_Kt6p",
        "TPS2_Kpi": None,
        "NTH1_Ktre": "nth_Ktrh",
        "HXK1_vmax": "hxk_Vmax",
        "PGI1_vmax": "pgi_Vmax",
        "PGM1_vmax": "pgm_Vmax",
        "TPS1_vmax": "tps_Vmax",
        "TPS2_vmax": "tpp_Vmax",
        "NTH1_vmax": "nth_Vmax",
    }
    # print(list(new_params.keys()))

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
    for para in fit_params:
        p, s, f, _t = optimize(
            smallbone,
            vhd_dataset,
            para,
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
        )
        save_traj(p, f"p_{para}_{time_stamp}.nc")
        save_traj(s, f"s_{para}_{time_stamp}.nc")
        save_traj(f, f"f_{para}_{time_stamp}.nc")


if __name__ == "__main__":
    main()
