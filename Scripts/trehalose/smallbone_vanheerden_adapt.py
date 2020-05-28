# -*- coding: utf-8 -*-
"""
actually run the Smallbone model on van Heerden's dataset
"""
# %%
from pyADAPT.examples import Smallbone2011
from pyADAPT.optimize import Optimizer, ITER, optimize
from van_heerden_preprocess import vhd_dataset
import pyADAPT.trajectory as traj
import numpy as np
import platform
from mat4py import loadmat
import os
import time


def main():
    matpath = os.path.join(
        os.path.dirname(__file__),
        os.pardir,
        os.pardir,
        r"data/trehalose/parametersDavid.mat",
    )
    time_stamp = time.strftime("%Y-%m-%d_%H.%M.%S", time.localtime())
    new_params = loadmat(matpath)["p"]
    syn = {
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
    print(list(new_params.keys()))

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

    if "compute" in platform.node():
        p, s, f, t = optimize(
            smallbone,
            vhd_dataset,
            *fit_params,
            #                  alternative initial paramters from david's email (not including udg)
            #                            (Glt)                      (TPS2)   (TPS1)          (!missing)
            #                  pgi_Vmax hxt_Vmax hxk_Vmax pgm_Vmax tpp_Vmax tps_Vmax nth_Vmax ugp_Vmax
            initial_parameters=[13.4667, 3.67, 4.75, 100, 81.45, 1000, 100, 36.8200],
            n_core=30,
            n_iter=120,
            delta_t=2,
            #               g1p g6p trh t6p udg
            weights=np.array([1, 1, 0.5, 1, 0.5]),
            timeout=100,
            attempt_limit=100,
            lambda_r=10,
            odesolver="Radau",
        )
        traj.save(p, f"p-{time_stamp}.nc")
        traj.save(s, f"s-{time_stamp}.nc")
        traj.save(f, f"f-{time_stamp}.nc")
    else:
        optim = Optimizer(
            model=smallbone, dataset=vhd_dataset, parameter_names=fit_params
        )
        print(optim.parameters)
        initial_parameters = [13.4667, 3.67, 4.75, 100, 81.45, 1000, 100, 36.8200]
        optim.parameters.loc[:, "init"] = initial_parameters

        # optim.run(n_core=1, n_iter=1, delta_t=2)
        print(optim.parameters)
        print(optim.state_mask)
        print(optim.model.state_order)
        print(optim.model.flux_order)


if __name__ == "__main__":
    main()
