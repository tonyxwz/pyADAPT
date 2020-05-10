"""
Steady state simulation of Smallbone model -- Classical approach
"""
#%%
from pyADAPT.io import read_data_raw, read_mat
from mat4py import loadmat
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from pyADAPT.dataset import DataSet, plot_splines

#%%
matpath = r"data\trehalose\vHeerden_trehalose_data_micromolgdw.mat"
vanHeerden = loadmat(matpath)

legenda_meta = np.squeeze(np.array(vanHeerden["data"]["legenda_metabolites"]))
meta = np.squeeze(np.array(vanHeerden["data"]["metabolites"]))
time_meta = np.squeeze(np.array(vanHeerden["data"]["time_metabolites"]))

legenda_nucleotides = np.squeeze(np.array(vanHeerden["data"]["legenda_nucleotides"]))
nucleotides = np.squeeze(np.array(vanHeerden["data"]["nucleotides"]))
time_nucleotides = np.squeeze(np.array(vanHeerden["data"]["time_nucleotides"]))

legenda_fluxes = np.squeeze(np.array(vanHeerden["data"]["legenda_fluxes"]))
fluxes = np.squeeze(np.array(vanHeerden["data"]["fluxes"]))
time_fluxes = np.squeeze(np.array(vanHeerden["data"]["time_fluxes"]))

totP = np.squeeze(np.array(vanHeerden["data"]["totP"]))
time_totP = np.squeeze(np.array(vanHeerden["data"]["time_totP"]))

# only use states, drop the fluxes

meta_std = meta * np.random.random_sample(meta.shape) * 0.2
# print(len(legenda_meta))
# print(legenda_meta)
# print(time_meta)
# print(len(meta))
# print(len(meta_std))

# create fake standard deviations
raw_meta = dict()
meta_specs = {
    "name": "van Heerden's dataset for trehalose cycle",
    "description": "For testing with simulation of Smallbone model",
    "raw_path": None,
    "groups": None,
}

meta_struct = dict()
raw_meta["t"] = time_meta

rename_map = {
    "G1P": "g1p",
    "G6P": "g6p",
    "Tre": "trh",
    "Tre6P": "t6p",
    "UDP_Glc": "udg"
}

for i, s in enumerate(legenda_meta):
    if not s in ['G1P', 'G6P', 'Tre', 'Tre6P', 'UDP_Glc']:
        continue
    s = rename_map[s]
    raw_meta[s + "_stds"] = meta_std[:, i]
    raw_meta[s + "_means"] = meta[:, i]
    tmp = dict()
    tmp["means"] = s + "_means"
    tmp["stds"] = s + "_stds"
    tmp["time"] = "t"
    tmp["unit"] = "mL"
    tmp["time_unit"] = "TBD"
    meta_struct[s] = tmp
meta_specs["structure"] = meta_struct
raw_meta["vhd"] = raw_meta
vhd_dataset = DataSet(name="vhd", raw_data=raw_meta, data_specs=meta_specs)

# fig, axes = plt.subplots(2, 3, figsize=(15, 8))
# plot_splines(
#     vhd_dataset, 100, 100, axes=axes
# )
# fig.tight_layout()
# # plt.show()

# We have the states now, and we can do a classical parameter estimation


# %%
from pyADAPT.examples import Smallbone2011
smallbone = Smallbone2011()

from pyADAPT.optimize import Optimizer, ITER, optimize

n_iter = 1
n_ts = 50

fit_params = list()

for id in smallbone.parameters.index:
    if id[-5:] == "_Vmax":
        fit_params.append(id)
fit_params

# optim = Optimizer(model=smallbone, dataset=vhd_dataset, parameter_names=fit_params)
# optim.run(n_core=1, n_iter=1, delta_t=2)
# %%
adapt_result = optimize(smallbone, vhd_dataset, *fit_params, n_core=1, n_iter=1, delta_t=2)