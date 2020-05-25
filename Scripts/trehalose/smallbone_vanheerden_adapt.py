# -*- coding: utf-8 -*-
"""
actually run the Smallbone model on van Heerden's dataset
"""
# %%
from pyADAPT.examples import Smallbone2011
from pyADAPT.optimize import Optimizer, ITER, optimize
from van_heerden_preprocess import vhd_dataset
import pyADAPT.trajectory as traj

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

    # %%

    p, s, f, t = optimize(
        smallbone, vhd_dataset, *fit_params, n_core=30, n_iter=120, delta_t=2, verbose=ITER
    )

    traj.save(p, "p.nc")
    traj.save(s, "s.nc")
    traj.save(f, "f.nc")


if __name__ == "__main__":
    main()