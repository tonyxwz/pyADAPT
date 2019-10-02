import matplotlib.pyplot as plt
import numpy as np

from pyADAPT import ADAPTapp, DataSet, Model
from toy_model import ToyModel


toy = ToyModel()
data = DataSet(raw_data_path='data/toyModel/toyData.npy',
        data_info_path='data/toyModel/toyData.yaml')
# print(isinstance(toy, Model))
app = ADAPTapp(toy, data)

# print(toy.constants)
# print(toy.parameters)
# print(toy.states)
# print(toy.observables)

# app.set_options(dict())
app.run(n_iter=2)
# app.set_options(hh=200)

x0 = toy.states
tspan = [-100, 0]
t_eval = np.linspace(*tspan, 100)

y = toy.compute_states(tspan, x0, t_eval=t_eval)

import matplotlib.pyplot as plt

plt.figure(1)
plt.plot(t_eval, y[0,:])
plt.plot(t_eval, y[1,:])
plt.plot(t_eval, y[2,:])
plt.plot(t_eval, y[3,:])

plt.legend(['s1', 's2', 's3', 's4'])
plt.figure(2)
plt.plot(t_eval, y[1,:])
plt.show()
