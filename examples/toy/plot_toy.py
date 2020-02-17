#%%
import matplotlib.pyplot as plt
import numpy as np

traj = np.load(
    r"C:\Users\tonyx\source\pyADAPT\examples\toy\traj-2020-02-16_15-59-21-046339.npy"
)
states = np.load(
    r"C:\Users\tonyx\source\pyADAPT\examples\toy\states-2020-02-16_15-59-21-046339.npy"
)

#%% Plot the trajectories

n_traj = traj.shape[0]
n_params = traj.shape[1]

n_plot = states.shape[0]
n_states = states.shape[1]


fig3, axes3 = plt.subplots(ncols=2, nrows=2)
for i in range(n_states):
    ax = axes3[i // 2, i % 2]
    for l in range(n_plot):
        ax.plot(states[l, i, 1:])
    ax.set_title("s" + str(i + 1))
fig3.tight_layout()

fig4, axes3 = plt.subplots()
for i in range(n_traj):
    axes3.plot(traj[i, 0, 1:])
axes3.set_title("k1")
fig4.tight_layout()

# %%
print(traj[:, :, 0])

plt.show()
# %%
