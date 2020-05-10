from pyADAPT.basemodel import BaseModel
import numpy as np


class Smallbone2011(BaseModel):
    def __init__(self):
        self.name = "Smallbone_Trehalose_Model"
        self.description = "Trehalose cycle"
        self.add_parameter("cell", 1.00000, False, 0.0, np.inf)
        self.add_parameter("medium", 1.00000, False, 0.0, np.inf)
        self.add_parameter("heat", 0.00000, False, 0.0, np.inf)
        self.add_parameter("glc_0", 0.09765, False, 0.0, np.inf)
        self.add_parameter("g1p_0", 0.10000, False, 0.0, np.inf)
        self.add_parameter("g6p_0", 2.67500, False, 0.0, np.inf)
        self.add_parameter("trh_0", 0.05000, False, 0.0, np.inf)
        self.add_parameter("t6p_0", 0.02000, False, 0.0, np.inf)
        self.add_parameter("udg_0", 0.70000, False, 0.0, np.inf)
        self.add_parameter("pgi_Vmax", 1071.00000, False, -np.inf, np.inf)
        self.add_parameter("pgi_Kg6p", 1.40000, False, -np.inf, np.inf)
        self.add_parameter("pgi_Kf6p", 0.29000, False, -np.inf, np.inf)
        self.add_parameter("pgi_Keq", 0.30000, False, -np.inf, np.inf)
        self.add_parameter("pgi_shock", 1.00000, False, -np.inf, np.inf)
        self.add_parameter("hxt_Vmax", 97.24000, False, -np.inf, np.inf)
        self.add_parameter("hxt_Kglc", 1.19180, False, -np.inf, np.inf)
        self.add_parameter("hxt_Ki", 0.91000, False, -np.inf, np.inf)
        self.add_parameter("hxt_shock", 8.00000, False, -np.inf, np.inf)
        self.add_parameter("hxk_Vmax", 289.60000, False, -np.inf, np.inf)
        self.add_parameter("hxk_Kglc", 0.08000, False, -np.inf, np.inf)
        self.add_parameter("hxk_Katp", 0.15000, False, -np.inf, np.inf)
        self.add_parameter("hxk_Kg6p", 30.00000, False, -np.inf, np.inf)
        self.add_parameter("hxk_Kadp", 0.23000, False, -np.inf, np.inf)
        self.add_parameter("hxk_Keq", 2000.00000, False, -np.inf, np.inf)
        self.add_parameter("hxk_Kit6p", 0.04000, False, -np.inf, np.inf)
        self.add_parameter("hxk_shock", 8.00000, False, -np.inf, np.inf)
        self.add_parameter("pgm_Vmax", 0.35450, False, -np.inf, np.inf)
        self.add_parameter("pgm_Kg6p", 0.05000, False, -np.inf, np.inf)
        self.add_parameter("pgm_Kg1p", 0.02300, False, -np.inf, np.inf)
        self.add_parameter("pgm_Keq", 0.16670, False, -np.inf, np.inf)
        self.add_parameter("pgm_shock", 16.00000, False, -np.inf, np.inf)
        self.add_parameter("tpp_Vmax", 6.50000, False, -np.inf, np.inf)
        self.add_parameter("tpp_Kt6p", 0.50000, False, -np.inf, np.inf)
        self.add_parameter("tpp_shock", 18.00000, False, -np.inf, np.inf)
        self.add_parameter("tps_Vmax", 1.37100, False, -np.inf, np.inf)
        self.add_parameter("tps_Kg6p", 3.80000, False, -np.inf, np.inf)
        self.add_parameter("tps_Kudg", 0.88600, False, -np.inf, np.inf)
        self.add_parameter("tps_shock", 12.00000, False, -np.inf, np.inf)
        self.add_parameter("tps_activity", 1.00000, False, -np.inf, np.inf)
        self.add_parameter("nth_Vmax", 15.20000, False, -np.inf, np.inf)
        self.add_parameter("nth_Ktrh", 2.99000, False, -np.inf, np.inf)
        self.add_parameter("nth_shock", 6.00000, False, -np.inf, np.inf)
        self.add_parameter("ugp_Vmax", 36.82000, False, -np.inf, np.inf)
        self.add_parameter("ugp_Kutp", 0.11000, False, -np.inf, np.inf)
        self.add_parameter("ugp_Kiutp", 0.11000, False, -np.inf, np.inf)
        self.add_parameter("ugp_Kg1p", 0.32000, False, -np.inf, np.inf)
        self.add_parameter("ugp_Kiudg", 0.00350, False, -np.inf, np.inf)
        self.add_parameter("ugp_shock", 16.00000, False, -np.inf, np.inf)

        # boundary conditions
        self.add_parameter('adp', 1.282000, False)
        self.add_parameter('atp', 2.525000, False)
        self.add_parameter('ppi', 1.000000, False)
        self.add_parameter('f6p', 0.625000, False)
        self.add_parameter('h', 1.000000, False)
        self.add_parameter('pho', 1.000000, False)
        self.add_parameter('udp', 0.281500, False)
        self.add_parameter('utp', 0.649100, False)
        self.add_parameter('h2o', 1.000000, False)
        self.add_parameter('glx', 100.000000, False)

        # self.add_state("glc", 0.09765, True)
        # self.add_state('g1p', 0.10000, True)
        # self.add_state('g6p', 2.67500, True)
        # self.add_state('trh', 0.05000, True)
        # self.add_state('t6p', 0.02000, True)
        # self.add_state('udg', 0.70000, True)

        self.sm = np.array([[0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
                            [-1.0, 0.0, 1.0, -1.0, 0.0, -1.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0]])

        super().__init__(
            state_order=['glc', 'g1p', 'g6p', 'trh', 't6p', 'udg'],
            initial_states=[0.09675, 0.1, 2.675, 0.05, 0.02, 0.7],
            flux_order=['pgi', 'hxt','hxk','pgm','tpp','tps','nth','ugp'],
        )

    def fluxes(self, t, x, p):
        # measured
        # wish:
        # glc = x['glc']
        glc = x[self['glc']]
        g1p = x[self['g1p']]
        g6p = x[self['g6p']]
        trh = x[self['trh']]
        t6p = x[self['t6p']]
        udg = x[self['udg']]

        # unmeasured and they don't change at all
        adp = p['adp']
        atp = p['atp']
        ppi = p['ppi']
        f6p = p['f6p']
        h = p['h']
        pho = p['pho']
        udp = p['udp']
        utp = p['utp']
        h2o = p['h2o']
        glx = p['glx']

        cell = p['cell']
        medium = p['medium']
        heat = p['heat']

        glc_0 = p['glc_0']
        g1p_0 = p['g1p_0']
        g6p_0 = p['g6p_0']
        trh_0 = p['trh_0']
        t6p_0 = p['t6p_0']
        udg_0 = p['udg_0']

        pgi_Vmax = p['pgi_Vmax']
        pgi_Kg6p = p['pgi_Kg6p']
        pgi_Kf6p = p['pgi_Kf6p']
        pgi_Keq = p['pgi_Keq']
        pgi_shock = p['pgi_shock']

        hxt_Vmax = p['hxt_Vmax']
        hxt_Kglc = p['hxt_Kglc']
        hxt_Ki = p['hxt_Ki']
        hxt_shock = p['hxt_shock']

        hxk_Vmax = p['hxk_Vmax']
        hxk_Kglc = p['hxk_Kglc']
        hxk_Katp = p['hxk_Katp']
        hxk_Kg6p = p['hxk_Kg6p']
        hxk_Kadp = p['hxk_Kadp']
        hxk_Keq = p['hxk_Keq']
        hxk_Kit6p = p['hxk_Kit6p']
        hxk_shock = p['hxk_shock']

        pgm_Vmax = p['pgm_Vmax']
        pgm_Kg6p = p['pgm_Kg6p']
        pgm_Kg1p = p['pgm_Kg1p']
        pgm_Keq = p['pgm_Keq']
        pgm_shock = p['pgm_shock']

        tpp_Vmax = p['tpp_Vmax']
        tpp_Kt6p = p['tpp_Kt6p']
        tpp_shock = p['tpp_shock']

        tps_Vmax = p['tps_Vmax']
        tps_Kg6p = p['tps_Kg6p']
        tps_Kudg = p['tps_Kudg']
        tps_shock = p['tps_shock']
        tps_activity = p['tps_activity']

        nth_Vmax = p['nth_Vmax']
        nth_Ktrh = p['nth_Ktrh']
        nth_shock = p['nth_shock']
        ugp_Vmax = p['ugp_Vmax']
        ugp_Kutp = p['ugp_Kutp']
        ugp_Kiutp = p['ugp_Kiutp']
        ugp_Kg1p = p['ugp_Kg1p']
        ugp_Kiudg = p['ugp_Kiudg']
        ugp_shock = p['ugp_shock']

        # G6P isomerase
        pgi = cell * pow(pgi_shock, heat) * pgi_Vmax / pgi_Kg6p * (
            g6p - f6p / pgi_Keq) / (1 + g6p / pgi_Kg6p + f6p / pgi_Kf6p)

        # glucose transport
        hxt = cell * pow(hxt_shock, heat) * hxt_Vmax * (
            glx - glc) / hxt_Kglc / (1 + (glx + glc) / hxt_Kglc +
                                     hxt_Ki * glx * glc / pow(hxt_Kglc, 2))
        # Hexokinase
        hxk = cell * pow(hxk_shock, heat) * hxk_Vmax / (
            hxk_Kglc * hxk_Katp) * (glc * atp - g6p * adp / hxk_Keq) / (
                (1 + glc / hxk_Kglc + g6p / hxk_Kg6p + t6p / hxk_Kit6p) *
                (1 + atp / hxk_Katp + adp / hxk_Kadp))
        # Phosphoglucomutase
        pgm = cell * pow(pgm_shock, heat) * pgm_Vmax / pgm_Kg6p * (
            g6p - g1p / pgm_Keq) / (1 + g6p / pgm_Kg6p + g1p / pgm_Kg1p)
        # T6P phosphatase
        tpp = cell * pow(
            tpp_shock, heat) * tpp_Vmax * t6p / tpp_Kt6p / (1 + t6p / tpp_Kt6p)
        # T6P synthase
        tps = cell * tps_activity * pow(
            tps_shock, heat) * tps_Vmax * g6p * udg / (tps_Kg6p * tps_Kudg) / (
                (1 + g6p / tps_Kg6p) * (1 + udg / tps_Kudg))
        # Trehalase
        nth = cell * pow(
            nth_shock, heat) * nth_Vmax * trh / nth_Ktrh / (1 + trh / nth_Ktrh)
        # UDP–glucose phosphorylase
        ugp = cell * pow(
            ugp_shock, heat) * ugp_Vmax * utp * g1p / (ugp_Kutp * ugp_Kg1p) / (
                ugp_Kiutp / ugp_Kutp + utp / ugp_Kutp + g1p / ugp_Kg1p +
                utp * g1p / (ugp_Kutp * ugp_Kg1p) +
                ugp_Kiutp / ugp_Kutp * udg / ugp_Kiudg + g1p * udg /
                (ugp_Kg1p * ugp_Kiudg))

        return np.array([pgi, hxt, hxk, pgm, tpp, tps, nth, ugp])

    def state_ode(self, t, x, p):
        # x: np.ndarray, how to use 'id' to index?
        v = self.fluxes(t, x, p)
        return np.dot(self.sm, v)

    def inputs(self, t):
        return super().inputs(t)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pyADAPT.dataset import DataSet
    from pyADAPT.optimize import optimize

    smallbone = Smallbone2011()
    # x0 = [
    #     0.09765, 0.10000, 2.67500, 0.05000, 0.02000, 0.70000, 1.28200, 2.52500,
    #     1.00000, 0.62500, 1.00000, 1.00000, 0.28150, 0.64910, 1.00000,
    #     100.00000
    # ]
    # t = np.linspace(0, 10, 1000)
    # x0[0] = 100  # glc
    # x0[-1] = 0.1  # glx
    # y = smallbone.compute_states(time_points=t, x0=x0)

    # y = y[[0, 1, 2, 3, 4, 5]]
    # for i in range(y.shape[0]):
    #     plt.plot(t, y[i, :])
    # plt.legend(['glc', 'g1p', 'g6p', 'trh', 't6p', 'udg'])
    # plt.title("glx=0.1")
    # plt.show()
    data = DataSet(raw_data_path="data/trehalose/smallbone2011_data.mat",
                   data_specs_path="data/trehalose/smallbone2011_data.yaml")
    ptraj, straj, time = optimize(smallbone,
                                  data,
                                  "tpp_Kt6p",
                                  n_iter=4,
                                  n_tstep=50,
                                  n_core=1,
                                  init_method=None,
                                  verbose=2)