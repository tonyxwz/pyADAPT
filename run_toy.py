import matplotlib.pyplot as plt
import numpy as np

from pyADAPT import ADAPT, DataSet, Model
from pyADAPT.models import ToyModel


if __name__ == "__main__":
    toy = ToyModel()
    data = DataSet(raw_data_path='data/toyModel/toyData.mat',
                data_specs_path='data/toyModel/toyData.yaml')
    # print(isinstance(toy, Model))
    adapt = ADAPT(toy, data)
    adapt.set_options(n_iter=20, n_ts=100)
    adapt.run(n_core=1)

    n_traj = adapt.trajectories.shape[0]
    n_params = adapt.trajectories.shape[1]
    fig, axes = plt.subplots(ncols=5, nrows=4, figsize=(15, 10))

    for i in range(n_traj):
        ax: plt.Axes = axes[i // 5, i % 5]
        keys = adapt.model.parameters.keys()

        for p in range(n_params):
            ax.plot(adapt.trajectories[i, p, :])
        ax.legend(keys)

    fig.tight_layout()

    fig2, axes2 = plt.subplots(ncols=5, nrows=4, figsize=(15, 10))
    n_plot = adapt.states.shape[0]
    n_states = adapt.states.shape[1]

    for i in range(n_plot):
        ax = axes2[i // 5, i % 5]
        for s in range(n_states):
            ax.plot(adapt.states[i, s, :])
        ax.legend(data.ordered_names)
    fig2.tight_layout()

    plt.show()

