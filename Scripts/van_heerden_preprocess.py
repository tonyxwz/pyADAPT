"""
return preprocessed van heerden dataset
1. use only states
2. divide by 2
3. pickout only ids which also appears in smallbone model
4. generate random standard deviation
5. pack into pyADAPT.DataSet

Steady state simulation of Smallbone model -- Classical approach

Metabolites
| state                 | id in Smallbone | id in vHeerden |
|-----------------------|-----------------|----------------|
| glucose 1-phosphate   | 1 g1p           | 12 G1P         |
| glucose 6-phosphate   | 2 g6p           | 1 G6P          |
| trehalose             | 3 trh           | 15 Tre         |
| trehalose 6-phosphate | 4 t6p           | 14 Tre6P       |
| UDP glucose           | 5 udg           | 13 UDP_Glc     |
Fluxes
| reaction                  | id in Smallbone | id in vHeerden         |
|---------------------------|-----------------|------------------------|
| T6P phosphatase           | 4 tpp           | (9) TPS2               |
| T6P synthase              | 5 tps           | (8) TPS1               |
| Trehalase                 | 6 nth           | (10) TREH              |
| UDP–glucose phosphorylase | 7 ugp           | (7) UPG                |
"""
#%%
from pyADAPT.io import read_data_raw, read_mat
from mat4py import loadmat
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from pyADAPT.dataset import DataSet, plot_splines

__all__ = ["vhd_dataset"]


def moving_average(data_set, periods=3):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, mode="valid")


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

# meta_std = meta * np.random.random_sample(meta.shape) * 0.05
# print(len(legenda_meta))
# print(legenda_meta)
# print(time_meta)
# print(len(meta))
# print(len(meta_std))

meta_std = meta * 0.01
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
    "UDP_Glc": "udg",
}

for i, s in enumerate(legenda_meta):
    if not s in ["G1P", "G6P", "Tre", "Tre6P", "UDP_Glc"]:
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


if __name__ == "__main__":
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    plot_splines(vhd_dataset, 100, 100, axes=axes)
    fig.tight_layout()
    plt.show()
