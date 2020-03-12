import matplotlib
import matplotlib.pyplot as plt


def create_axes(axes):
    if axes is None:
        fig, axes = plt.subplots()
    else:
        fig = None

    return fig, axes


def plot_error(state, axes=None):
    # TODO
    fig, axes = create_axes(axes)

    if fig is not None:
        plt.show()