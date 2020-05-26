# plot the success iterations of the simulation

#%% success iters
import numpy as np

i = np.loadtxt("success-iters")
print(i)
# %%
import pyADAPT.trajectory as traj

s = traj.load("s_5-22_smallbone.nc")
f = traj.load("f_5-22_smallbone.nc")
p = traj.load("p_5-22_smallbone.nc")

# %%
s = s.sel(iter=i)
f = f.sel(iter=i)
p = p.sel(iter=i)

# %%
import matplotlib.pyplot as plt

plt.style.use(["science", "grid"])
fig, axes = plt.subplots(2, 4, figsize=(9, 5))
traj.plot(p, axes=axes, color="blue", alpha=0.2)
traj.plot_mean(p, axes=axes, color="red")
fig.tight_layout()
fig.savefig("p_success.png", dpi=200)

# %%
# fig2 = plt.figure(figsize=(9, 5))
# ax1 = plt.subplot2grid((2, 6), loc=(0, 0), colspan=2)
# ax2 = plt.subplot2grid((2, 6), (0, 2), colspan=2)
# ax3 = plt.subplot2grid((2, 6), (0, 4), colspan=2)
# ax4 = plt.subplot2grid((2, 6), (1, 1), colspan=2)
# ax5 = plt.subplot2grid((2, 6), (1, 3), colspan=2)
# axes2 = np.asarray([ax1, ax2, ax3, ax4, ax5])
fig2, axes2 = plt.subplots(2, 3, figsize=(9, 5))
traj.plot(s, axes=axes2, color="blue", alpha=0.2)
traj.plot_mean(s, axes=axes2, color="red")
fig2.tight_layout()
fig2.savefig("s_success.png", dpi=200)


# %%
fig3, axes3 = plt.subplots(2, 4, figsize=(9, 5))
traj.plot(f, axes=axes3, color="blue", alpha=0.2)
traj.plot_mean(f, axes=axes3, color="red")
fig3.tight_layout()
fig3.savefig("f_success.png", dpi=200)


# %%
