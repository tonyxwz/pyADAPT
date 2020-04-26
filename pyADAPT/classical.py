"""
classic approach parameter estimation

mismatch between the different metabolites in Smallbone and vanHeerden
This metabolites have little influence on the model (as far as we know). That is why the experimentalis did not take the work to measure. I would fix the initial contentration to the same metabolite value the Smallbone 2011 model uses, and let it run to steady state.

Then, when you do least squares optimization, you can consider if you want to include it in the cost function or not. At least in the matlab ADAPT version it can be selected if a metabolite is observable or not.

In a way, you tell the algorithm to model it, starting at the literature value, but then to forget about it for error minimization.

1. initialize model
2. assign literature parameter
3. simulate until steady states on the time range of van Heerden's data set
4.
"""
import numpy as np
import matplotlib.pyplot as plt


def cpe(model, data, params):
    """ classical parameter estimation

    """


if __name__ == "__main__":
    from pyADAPT.examples import Smallbone2011

    model = Smallbone2011()

    time_points = np.linspace(0, 340)  # van heerden time range from 0 to 340
    x0 = [.09675, .1, 2.675, .05, .02, .7]  # Smallbone literature
    states = model.compute_states(time_points=time_points, x0=x0)

    plt.plot(time_points, states.T)
    plt.title('test')
    plt.show()
