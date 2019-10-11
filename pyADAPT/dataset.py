import numpy as np

from pyADAPT.core import State
from pyADAPT.io import read_data_info, read_data_raw


class DataSet(list):
    """
    dataset for ADAPT. containing phenotypes measured at different stages after
    intervention. list of States

    Job division from DataSet and Model: DataSet only stores the source data and
    yields new interpolant data. Model stores the interpolants as the
    trajectories.
    """
    def __init__(self, raw_data_path="", data_specs_path="",
                 raw_data={}, data_specs={}):
        """
        raw_data: phenotypes organized into a dictionary

        data_info: instructions of the data, such as which time variable 
            should be used for which state.
        """
        if data_specs:
            self.data_specs = data_specs
        else:
            self.data_specs = read_data_info(data_specs_path)

        if raw_data:
            self.raw_data = raw_data
        else:
            # group: to solve the stupid matlab issue
            group = self.data_specs['mat_group'] if 'mat_group' in self.data_specs else ""
            self.raw_data = read_data_raw(raw_data_path, group=group)

        self.structure = {}
        self.ordered_names = []

        for k,v in self.data_specs.items():
            self.__setattr__(k, v)

        for k, v in self.structure.items():
            # !print(v)
            time = self.raw_data[v['time']]
            means = self.raw_data[v['means']]
            stds = self.raw_data[v['stds']]
            s = State(name=k, time=time, means=means, stds=stds)
            # self[k] = s
            self.append(s)
            self.ordered_names.append(k)

    def interpolate(self, n_ts=100, method='Hermite'):
        """In every ADAPT iteration, this function is called once to get a new
        spline for the optimizer to fit (from t0 till the end). the length of
        the list of the splines should equal the number of states in the data.

        return
        ------
        numpy.ndarray in the same order as `self.ordered_names`
        ```
        [
            [
                [value, std],
                ......
                [value, std]
            ]
        ]
        ```
        """
        # TODO: add different interp methods
        # TODO: take care of states using different time points, t1, t2
        ainter_p = np.zeros((len(self), n_ts, 2))  #a stands for array
        # v = list(self.values())
        for i in range(len(self)):
            ainter_p[i, :, 0] = self[i].interp_values(n_ts=n_ts)
            ainter_p[i, :, 1] = self[i].interp_stds(n_ts=n_ts)
        return ainter_p

    def __getitem__(self, index):
        if type(index) is str:
            index = self.ordered_names.index(index)
        return super().__getitem__(index)


if __name__ == "__main__":

    from pprint import pprint, pformat
    import matplotlib.pyplot as plt
    # import seaborn as sns

    D = DataSet(raw_data_path='data/toyModel/toyData.npy',
        data_specs_path='data/toyModel/toyData.yaml')
    idp = D.interpolate(n_ts=10)

    # all for s1:
    pprint(idp[0, :, :])

    # select all the data from time step 3
    pprint(idp[:,3,:])

    # all for values
    pprint(idp[:, :, 0])

    n_interp = 100
    n_ts = 200
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # 4 interpolations, 4 states
    fig.canvas.set_window_title(f'DataSet Interpolation')
    fig.suptitle('Toy data interpolation (100)')
    np.random.seed(725)
    # plot all interpolations of values
    for i_interp in range(n_interp):
        idp = D.interpolate(n_ts=n_ts)
        for i_state in range(idp.shape[0]):
            ax:plt.Axes = axes[i_state//2, i_state%2]
            state:State = D[i_state]
            t_ = state.time
            t = np.linspace(t_[0], t_[-1], n_ts)
            ax.plot(t, idp[i_state, :, 0], color='red', alpha=0.15)  # values of s1
            ax.plot(t_, state.sampled_values, '.g', alpha=0.5, markersize=5)

            if not ax.title.get_label():
                ax.set_title(f"{D.ordered_names[i_state]}")
                ax.set_xlabel('days')
                err = ax.errorbar(t_, [d.mean for d in state],
                        yerr=[d.std for d in state],
                        fmt='.b',
                        uplims=True,
                        lolims=True)

    # fig.tight_layout()
    plt.show()
