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

    optim = Optimizer(model=smallbone, dataset=vhd_dataset, parameter_names=fit_params)

    # optim.run(n_core=1, n_iter=1, delta_t=2)
    print(optim.parameters)
    print(optim.state_mask)
    print(optim.model.state_order)
    print(optim.model.flux_order)

    if "compute" in platform.node():
        p, s, f, t = optimize(
            smallbone,
            vhd_dataset,
            *fit_params,
            #                  alternative initial paramters from david's email
            #
            initial_parameters=None,
            n_core=30,
            n_iter=120,
            delta_t=2,
            #               g1p g6p trh t6p udg
            weights=np.array([1, 1, 0.5, 1, 0.5]),
            timeout=100,
            attempt_limit=100
        )

        traj.save(p, "p.nc")
        traj.save(s, "s.nc")
        traj.save(f, "f.nc")


if __name__ == "__main__":
    main()
