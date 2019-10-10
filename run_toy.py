import matplotlib.pyplot as plt
import numpy as np

from pyADAPT import ADAPTapp, DataSet, Model
from toy_model import ToyModel
from logging import DEBUG

toy = ToyModel()
data = DataSet(raw_data_path='data/toyModel/toyData.npy', 
    data_info_path='data/toyModel/toyData.yaml')
# print(isinstance(toy, Model))
app = ADAPTapp(toy, data)
app.set_options(n_iter=2, n_ts=100, log_level=DEBUG)
# print(toy.constants)
# toy.parameters.pretty_print()
# print(toy.states)
# print(toy.observables)

app.randomize_parameters()
# toy.parameters.pretty_print()
# print([p.vary for p in toy.parameters.values()])
# print(toy.parameters.valuesdict())
app.run()

n_traj = len(app.trajectories)
axes:list
fig, axes = plt.subplots(ncols=3, nrows=1)

for i in range(n_traj):
    ax:plt.Axes = axes[i]
    keys = app.trajectories[i].keys()

    for v in app.trajectories[i].values():
        ax.plot(v)
    ax.legend(keys)
    
x0 = toy.states
tspan = [-100, 0]
t_eval = np.linspace(*tspan, 100)

plt.show()
# y = toy.compute_states(tspan, x0, t_eval=t_eval)

# import matplotlib.pyplot as plt

# plt.figure(1)
# plt.plot(t_eval, y[0,:])
# plt.plot(t_eval, y[1,:])
# plt.plot(t_eval, y[2,:])
# plt.plot(t_eval, y[3,:])

# plt.legend(['s1', 's2', 's3', 's4'])
# plt.figure(2)
# plt.plot(t_eval, y[1,:])
# plt.show()
