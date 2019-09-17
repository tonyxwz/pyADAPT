from pyADAPT import DataSet, Model, ADAPTapp

import matplotlib.pyplot as plt
import numpy as np

class ToyModel(Model):
    @staticmethod
    def odefunc(t, x, p, u):
        '''
        define odes as in natal's github repo. takes t: time, x: states,
        p: parameters and u: constant inputs as the arguments. return the
        differential of each state variable.
        '''
        u1 = u[0]
        u2 = u[1]

        s1 = x[0]
        s2 = x[1]
        s3 = x[2]
        s4 = x[3]

        k1 = p[0]
        k2 = p[1]
        k3 = p[2]
        k4 = p[3]
        k5 = p[4]

        v1 = k1 * u1 * s2
        v2 = k2 * u2 * s3
        v3 = k3 * s1
        v4 = k4 * s1
        v5 = k5 * s4
        ds3dt = v1 - v2

        dxdt = np.zeros(len(x))
        dxdt[0] = v1 - v3 - v4
        dxdt[1] = -v1 + v2
        dxdt[2] = ds3dt
        dxdt[3] = v4 - v5

        return dxdt

D = DataSet(raw_data_path='data/toyModel/toyData.npy',
    data_info_path='data/toyModel/toyData.yaml')
spl = D.generate_splines()
which_plot = 'std'
# plt.subplot(211)
plt.plot(spl['s1'][which_plot])
plt.plot(spl['s2'][which_plot])
plt.plot(spl['s3'][which_plot])
plt.plot(spl['s4'][which_plot])
# plt.legend(['s1', 's2', 's3', 's4'])

which_plot = 'value'
# plt.subplot(212)
plt.plot(spl['s1'][which_plot])
plt.plot(spl['s2'][which_plot])
plt.plot(spl['s3'][which_plot])
plt.plot(spl['s4'][which_plot])
plt.legend(['s1', 's2', 's3', 's4'] + [s+"_"+which_plot for s in ['s1', 's2', 's3', 's4']])
plt.show()