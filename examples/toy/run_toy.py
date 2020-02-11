import matplotlib.pyplot as plt
import numpy as np

from pyADAPT import ADAPT, DataSet, Model
from toy_model import ToyModel


if __name__ == "__main__":
    toy = ToyModel()
    # data = DataSet(raw_data_path=r"D:\Weizhou\ADAPT\AMF\models\toyModel\paperData.mat",
    data = DataSet(raw_data_path="data/toyModel/toyData.mat",
                data_specs_path='data/toyModel/toyData.yaml')
    # print(isinstance(toy, Model))
    adapt = ADAPT(toy, data)
    adapt.set_options(n_iter=20, n_ts=100, smin=-2, smax=2, ss_time=1000, seed=124)
    adapt.run(n_core=4)

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
            ax.plot(adapt.states[i, s, 1:])
        ax.legend(data.ordered_names)
    fig2.tight_layout()

    fig3, axes3 = plt.subplots(ncols=2, nrows=2)
    for i in range(n_states):
        ax = axes3[i//2, i%2]
        for l in range(n_plot):
            ax.plot(adapt.states[l, i, 1:])
        ax.set_title('s'+str(i+1))
    fig3.tight_layout()

    fig4, axes3 = plt.subplots()
    for i in range(n_traj):
        axes3.plot(adapt.trajectories[i, 0, 1:])
    axes3.set_title('k1')
    fig4.tight_layout()

    plt.show()
