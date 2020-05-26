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


def main():
    smallbone = Smallbone2011()

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
            initial_parameters=[13.4667,    3.67,    4.75,     100,   81.45,   1000,     100,  36.8200],
            n_core=30,
            n_iter=120,
            delta_t=2,
            #               g1p g6p trh t6p udg
            weights=np.array([1, 1, 0.5, 1, 0.5]),
            timeout=100,
            attempt_limit=100,
            lambda_r=10
        )
        traj.save(p, "p.nc")
        traj.save(s, "s.nc")
        traj.save(f, "f.nc")
    else:
        optim = Optimizer(model=smallbone, dataset=vhd_dataset, parameter_names=fit_params)
        print(optim.parameters)
        initial_parameters=[13.4667,    3.67,    4.75,     100,   81.45,   1000,     100,  36.8200]
        optim.parameters.loc[:, 'init'] = initial_parameters

        # optim.run(n_core=1, n_iter=1, delta_t=2)
        print(optim.parameters)
        print(optim.state_mask)
        print(optim.model.state_order)
        print(optim.model.flux_order)

if __name__ == "__main__":
    main()
