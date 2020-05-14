import pyADAPT.trajectory as traj
import matplotlib.pyplot as plt
import numpy as np


states = traj.load("s.nc")
params = traj.load("p.nc")
fluxes = traj.load("f.nc")
print(states.shape)
print(params.shape)
print(fluxes.shape)
# fig, axes = plt.subplots(2,3)
# traj.plot(states, axes=axes, color="red", alpha=0.1)
# plt.show()

fig, axes = plt.subplots()

# axes.plot(states[0, :, 0])
correfs = np.zeros(120)
for iter in range(120):
    correfs[iter] = np.correlate(states[0, :, 0], states[iter, :, 0])
axes.plot(correfs)
axes.set_title('Correlation')

fig2, axes2 = plt.subplots(2, 2)
for iter in range(120):
    axes2[(iter // 30) // 2, (iter // 30) % 2].plot(
        states.coords["time"], states[iter, :, 2]
    )
fig2.suptitle(states.coords['id'][2].data[()])
# plt.show()

fig.savefig("fig1.png")
fig2.savefig("fig2.png")
