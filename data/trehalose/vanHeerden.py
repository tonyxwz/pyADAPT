# script to create an ADAPT dataset from van Heerden's data
import numpy as np
from mat4py import loadmat
vanHeerden_mat = loadmat(
    'data/trehalose/vHeerden_trehalose_data_micromolgdw.mat')

legenda_meta = np.squeeze(
    np.array(vanHeerden_mat['data']["legenda_metabolites"]))
meta = np.squeeze(np.array(vanHeerden_mat['data']["metabolites"]))
time_meta = np.squeeze(np.array(vanHeerden_mat['data']["time_metabolites"]))

legenda_nucleotides = np.squeeze(
    np.array(vanHeerden_mat['data']["legenda_nucleotides"]))
nucleotides = np.squeeze(np.array(vanHeerden_mat['data']["nucleotides"]))
time_nucleotides = np.squeeze(
    np.array(vanHeerden_mat['data']["time_nucleotides"]))

legenda_fluxes = np.squeeze(np.array(vanHeerden_mat['data']["legenda_fluxes"]))
fluxes = np.squeeze(np.array(vanHeerden_mat['data']["fluxes"]))
time_fluxes = np.squeeze(np.array(vanHeerden_mat['data']["time_fluxes"]))

# TODO: ask david what is totP?
totP = np.squeeze(np.array(vanHeerden_mat['data']['totP']))
time_totP = np.squeeze(np.array(vanHeerden_mat['data']['time_totP']))

# Metabolites
# | state                 | id in Smallbone | id in vHeerden |
# |-----------------------|-----------------|----------------|
# | glucose 1-phosphate   | 1 g1p           | 12 G1P         |
# | glucose 6-phosphate   | 2 g6p           | 1 G6P          |
# | trehalose             | 3 trh           | 15 Tre         |
# | trehalose 6-phosphate | 4 t6p           | 14 Tre6P       |
# | UDP glucose           | 5 udg           | 13 UDP_Glc     |
# Fluxes
# | reaction                  | id in Smallbone | id in vHeerden         |
# |---------------------------|-----------------|------------------------|
# | T6P phosphatase           | 4 tpp           | (9) TPS2               |
# | T6P synthase              | 5 tps           | (8) TPS1               |
# | Trehalase                 | 6 nth           | (10) TREH              |
# | UDP–glucose phosphorylase | 7 ugp           | (7) UPG                |
