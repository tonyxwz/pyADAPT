import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numpy.random import rand, randn
from ode_true import ode_true
from fluxes_true import fluxes_true



plt.figure()
im = plt.imread('true_system.png')
plt.axis('off')
plt.imshow(im)
x0 = [1.03, .38, .62, .52, .52]
tspan = [-1e3, 0]

k6range = np.logspace(-2, -3.5, 5)
k60 = k6range[0]

ydata = []
v = []
for k6 in k6range:
    sol = solve_ivp(lambda t, x: ode_true(t, x, k6),
                tspan, x0)
    # print(sol.y[:,-1])
    ydata.append(sol.y[:,-1])
    v.append(fluxes_true(0, ydata[-1]))

ydata = np.array(ydata)
v = np.array(v)

plt.figure()
plt.plot(range(5), v, "*-")
plt.legend(["v1", "v2", "v3", "v4", "v5"])
plt.title('fluxes')

np.random.seed(100)
stddata = .2 * np.diag(ydata[0,:]).dot(rand(*ydata.shape))

ydataNoisy = ydata  + stddata * randn(*ydata.shape)

concdata = ydataNoisy[:, 0:4]
proteindata = np.abs(ydataNoisy[:, 4])
fluxes = v

fig, axes = plt.subplots(2, 2)

for i in range(4):
    ax = axes[i//2, i%2]
    ax.bar(np.arange(5), concdata[:,i])
    ax.set_title(f"$S_{i+1}$")
    ax.set_xticks(np.arange(5))
plt.tight_layout()


plt.figure()
plt.bar(np.arange(5), proteindata)
plt.title('$R_1$')


# short-term kinetic response

for k6 in k6range:
    tspan1 = [-1e3, 5]
    sol1 = solve_ivp(lambda t, x: ode_true(t, x, k6),
                     tspan1, x0)
    tspan2 = [5, 15];
    sol2 = solve_ivp(lambda t, x: ode_true(t, x, k6),
                     tspan2, sol1.y[:, -1])
    y = np.concatenate( (sol1.y[0:4, :], sol2.y[0:4, 1:]), axis=1)
    t = np.concatenate((sol1.t, sol2.t[1:]))
    plt.figure()
    plt.plot(t, y.T)
    plt.xlim(0, sol2.t[-1])
    plt.legend(['S1', 'S2', 'S3', 'S4'])
    plt.title(f"$k_6={k6}$")



