# plot the success iterations of the simulation

#%% success iters
import numpy as np

i = np.loadtxt("success-iters2")
print(i)
# %%
import pyADAPT.visualize as vis
from pyADAPT.io import load_traj, save_traj

s = load_traj("s_sbml-params_2020-05-28_19.44.30.nc")
f = load_traj("f_sbml-params_2020-05-28_19.44.30.nc")
p = load_traj("p_sbml-params_2020-05-28_19.44.30.nc")

# %%
s = s.sel(iter=i)
f = f.sel(iter=i)
p = p.sel(iter=i)

# %%
import matplotlib.pyplot as plt

plt.style.use(["science", "grid"])
fig, axes = plt.subplots(2, 4, figsize=(9, 5))
vis.plot(p, axes=axes, color="blue", alpha=0.2)
vis.plot_mean(p, axes=axes, color="red")
fig.tight_layout()
fig.savefig("p_success2.png", dpi=200)

# %%
# fig2 = plt.figure(figsize=(9, 5))
# ax1 = plt.subplot2grid((2, 6), loc=(0, 0), colspan=2)
# ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
# ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
# ax4 = plt.subplot2grid((2, 6), (1, 1), colspan=2)
# ax5 = plt.subplot2grid((2, 6), (1, 3), colspan=2)
# axes2 = np.asarray([ax1, ax2, ax3, ax4, ax5])
fig2, axes2 = plt.subplots(2, 3, figsize=(9, 5))
vis.plot(s, axes=axes2, color="blue", alpha=0.2)
vis.plot_mean(s, axes=axes2, color="red")
fig2.tight_layout()
fig2.savefig("s_success2.png", dpi=200)


# %%
fig3, axes3 = plt.subplots(2, 4, figsize=(9, 5))
vis.plot(f, axes=axes3, color="blue", alpha=0.2)
vis.plot_mean(f, axes=axes3, color="red")
fig3.tight_layout()
fig3.savefig("f_success2.png", dpi=200)


# %%
