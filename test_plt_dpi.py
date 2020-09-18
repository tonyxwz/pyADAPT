import matplotlib.pyplot as plt
import numpy as np

t = np.linspace(-np.pi, np.pi, 100)
plt.plot(t, np.tanh(t))
plt.show()

