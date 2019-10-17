import matplotlib.pyplot as plt
import numpy as np

from pyADAPT import ADAPT, DataSet, Model
from pyADAPT.models import ToyModel
from logging import DEBUG

toy = ToyModel()
data = DataSet(raw_data_path='data/toyModel/toyData.npy', 
    data_specs_path='data/toyModel/toyData.yaml')
# print(isinstance(toy, Model))
adapt = ADAPT(toy, data)
adapt.set_options(n_iter=20, n_ts=100, log_level=DEBUG)
# print(toy.constants)
# toy.parameters.pretty_print()
# print(toy.states)
# print(toy.observables)

# toy.parameters.pretty_print()
# print([p.vary for p in toy.parameters.values()])
# print(toy.parameters.valuesdict())
adapt.run()

n_traj = adapt.trajectories.shape[0]
n_params = adapt.trajectories.shape[1]
fig, axes = plt.subplots(ncols=5, nrows=4, figsize=(15,10))

for i in range(n_traj):
    
    ax:plt.Axes = axes[i//5, i%5]
    keys = adapt.model.parameters.keys()

    for p in range(n_params):
        ax.plot(adapt.trajectories[i, p, :])
    ax.legend(keys)

x0 = toy.states
tspan = [-100, 0]
t_eval = np.linspace(*tspan, 100)
fig.tight_layout()
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
