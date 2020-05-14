# Classical parameter estimation of Smallbone model
#   using the data from van Heerden's publication

from pyADAPT.io import read_data_raw, read_mat
from mat4py import loadmat
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

matpath = r'data\trehalose\vHeerden_trehalose_data_micromolgdw.mat'
vanHeerden = loadmat(matpath)

legenda_meta = np.squeeze(np.array(vanHeerden['data']["legenda_metabolites"]))
meta = np.squeeze(np.array(vanHeerden['data']["metabolites"]))
time_meta = np.squeeze(np.array(vanHeerden['data']["time_metabolites"]))

legenda_nucleotides = np.squeeze(
    np.array(vanHeerden['data']["legenda_nucleotides"]))
nucleotides = np.squeeze(np.array(vanHeerden['data']["nucleotides"]))
time_nucleotides = np.squeeze(np.array(vanHeerden['data']["time_nucleotides"]))

legenda_fluxes = np.squeeze(np.array(vanHeerden['data']["legenda_fluxes"]))
fluxes = np.squeeze(np.array(vanHeerden['data']["fluxes"]))
time_fluxes = np.squeeze(np.array(vanHeerden['data']["time_fluxes"]))

print(np.squeeze(np.array(meta)))

fig1, axes = plt.subplots()
for m in range(meta.shape[1]):
    axes.plot(time_meta, meta[:, m])
axes.legend(legenda_meta)
axes.set_title('meta')

fig2, axes = plt.subplots()
for m in range(nucleotides.shape[1]):
    axes.plot(time_nucleotides, nucleotides[:, m])
axes.legend(legenda_nucleotides)
axes.set_title('nucleotides')

fig3, axes = plt.subplots()
for m in range(fluxes.shape[1]):
    axes.plot(time_fluxes, fluxes[:, m])
axes.legend(legenda_fluxes)
axes.set_title('fluxes')

fig4, axes = plt.subplots()
axes.plot(vanHeerden['data']['time_totP'], vanHeerden["data"]['totP'])
axes.set_title('totP')

plt.show()
