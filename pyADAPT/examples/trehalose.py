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

        self.add_state("glc", -1954.423351, True)
        self.add_state('g1p', 0.058932, True)
        self.add_state('g6p', -8778.849594, True)
        self.add_state('trh', 0.071574, True)
        self.add_state('t6p', 0.028853, True)
        self.add_state('udg', 0.309360, True)
        self.add_state('adp', 1.282000, True)
        self.add_state('atp', 2.525000, True)
        self.add_state('ppi', 1.000000, True)
        self.add_state('f6p', 0.625000, True)
        self.add_state('h', 1.000000, True)
        self.add_state('pho', 1.000000, True)
        self.add_state('udp', 0.281500, True)
        self.add_state('utp', 0.649100, True)
        self.add_state('h2o', 1.000000, True)
        self.add_state('glx', 100.000000, True)

        self.sm = np.array([[0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
                            [-1.0, 0.0, 1.0, -1.0, 0.0, -1.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
        super().__init__()

    def fluxes(self, t, x, p):
        glc = x[0]
        g1p = x[1]
        g6p = x[2]
        trh = x[3]
        t6p = x[4]
        udg = x[5]
        adp = x[6]
        atp = x[7]
        ppi = x[8]
        f6p = x[9]
        h = x[10]
        pho = x[11]
        udp = x[12]
        utp = x[13]
        h2o = x[14]
        glx = x[15]

        cell = self.parameters.loc['cell', 'value']
        medium = self.parameters.loc['medium', 'value']
        heat = self.parameters.loc['heat', 'value']

        glc_0 = self.parameters.loc['glc_0', 'value']
        g1p_0 = self.parameters.loc['g1p_0', 'value']
        g6p_0 = self.parameters.loc['g6p_0', 'value']
        trh_0 = self.parameters.loc['trh_0', 'value']
        t6p_0 = self.parameters.loc['t6p_0', 'value']
        udg_0 = self.parameters.loc['udg_0', 'value']

        pgi_Vmax = self.parameters.loc['pgi_Vmax', 'value']
        pgi_Kg6p = self.parameters.loc['pgi_Kg6p', 'value']
        pgi_Kf6p = self.parameters.loc['pgi_Kf6p', 'value']
        pgi_Keq = self.parameters.loc['pgi_Keq', 'value']
        pgi_shock = self.parameters.loc['pgi_shock', 'value']

        hxt_Vmax = self.parameters.loc['hxt_Vmax', 'value']
        hxt_Kglc = self.parameters.loc['hxt_Kglc', 'value']
        hxt_Ki = self.parameters.loc['hxt_Ki', 'value']
        hxt_shock = self.parameters.loc['hxt_shock', 'value']

        hxk_Vmax = self.parameters.loc['hxk_Vmax', 'value']
        hxk_Kglc = self.parameters.loc['hxk_Kglc', 'value']
        hxk_Katp = self.parameters.loc['hxk_Katp', 'value']
        hxk_Kg6p = self.parameters.loc['hxk_Kg6p', 'value']
        hxk_Kadp = self.parameters.loc['hxk_Kadp', 'value']
        hxk_Keq = self.parameters.loc['hxk_Keq', 'value']
        hxk_Kit6p = self.parameters.loc['hxk_Kit6p', 'value']
        hxk_shock = self.parameters.loc['hxk_shock', 'value']

        pgm_Vmax = self.parameters.loc['pgm_Vmax', 'value']
        pgm_Kg6p = self.parameters.loc['pgm_Kg6p', 'value']
        pgm_Kg1p = self.parameters.loc['pgm_Kg1p', 'value']
        pgm_Keq = self.parameters.loc['pgm_Keq', 'value']
        pgm_shock = self.parameters.loc['pgm_shock', 'value']

        tpp_Vmax = self.parameters.loc['tpp_Vmax', 'value']
        tpp_Kt6p = self.parameters.loc['tpp_Kt6p', 'value']
        tpp_shock = self.parameters.loc['tpp_shock', 'value']

        tps_Vmax = self.parameters.loc['tps_Vmax', 'value']
        tps_Kg6p = self.parameters.loc['tps_Kg6p', 'value']
        tps_Kudg = self.parameters.loc['tps_Kudg', 'value']
        tps_shock = self.parameters.loc['tps_shock', 'value']
        tps_activity = self.parameters.loc['tps_activity', 'value']

        nth_Vmax = self.parameters.loc['nth_Vmax', 'value']
        nth_Ktrh = self.parameters.loc['nth_Ktrh', 'value']
        nth_shock = self.parameters.loc['nth_shock', 'value']
        ugp_Vmax = self.parameters.loc['ugp_Vmax', 'value']
        ugp_Kutp = self.parameters.loc['ugp_Kutp', 'value']
        ugp_Kiutp = self.parameters.loc['ugp_Kiutp', 'value']
        ugp_Kg1p = self.parameters.loc['ugp_Kg1p', 'value']
        ugp_Kiudg = self.parameters.loc['ugp_Kiudg', 'value']
        ugp_shock = self.parameters.loc['ugp_shock', 'value']

        pgi = cell * pow(pgi_shock, heat) * pgi_Vmax / pgi_Kg6p * (
            g6p - f6p / pgi_Keq) / (1 + g6p / pgi_Kg6p + f6p / pgi_Kf6p)
        hxt = cell * pow(hxt_shock, heat) * hxt_Vmax * (
            glx - glc) / hxt_Kglc / (1 + (glx + glc) / hxt_Kglc +
                                     hxt_Ki * glx * glc / pow(hxt_Kglc, 2))
        hxk = cell * pow(hxk_shock, heat) * hxk_Vmax / (
            hxk_Kglc * hxk_Katp) * (glc * atp - g6p * adp / hxk_Keq) / (
                (1 + glc / hxk_Kglc + g6p / hxk_Kg6p + t6p / hxk_Kit6p) *
                (1 + atp / hxk_Katp + adp / hxk_Kadp))
        pgm = cell * pow(pgm_shock, heat) * pgm_Vmax / pgm_Kg6p * (
            g6p - g1p / pgm_Keq) / (1 + g6p / pgm_Kg6p + g1p / pgm_Kg1p)
        tpp = cell * pow(
            tpp_shock, heat) * tpp_Vmax * t6p / tpp_Kt6p / (1 + t6p / tpp_Kt6p)
        tps = cell * tps_activity * pow(
            tps_shock, heat) * tps_Vmax * g6p * udg / (tps_Kg6p * tps_Kudg) / (
                (1 + g6p / tps_Kg6p) * (1 + udg / tps_Kudg))
        nth = cell * pow(
            nth_shock, heat) * nth_Vmax * trh / nth_Ktrh / (1 + trh / nth_Ktrh)
        ugp = cell * pow(
            ugp_shock, heat) * ugp_Vmax * utp * g1p / (ugp_Kutp * ugp_Kg1p) / (
                ugp_Kiutp / ugp_Kutp + utp / ugp_Kutp + g1p / ugp_Kg1p +
                utp * g1p / (ugp_Kutp * ugp_Kg1p) +
                ugp_Kiutp / ugp_Kutp * udg / ugp_Kiudg + g1p * udg /
                (ugp_Kg1p * ugp_Kiudg))

        return np.array([pgi, hxt, hxk, pgm, tpp, tps, nth, ugp])

    def odefunc(self, t, x, p):
        v = self.fluxes(t, x, p)
        self.flux_trajectory.append(v)
        return np.dot(self.sm, v)

    def inputs(self, t):
        return super().inputs(t)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    smallbone = Smallbone2011()
    x0 = [
        0.09765, 0.10000, 2.67500, 0.05000, 0.02000, 0.70000, 1.28200, 2.52500,
        1.00000, 0.62500, 1.00000, 1.00000, 0.28150, 0.64910, 1.00000,
        100.00000
    ]
    t = np.linspace(0, 10, 1000)
    x0[0] = 100  # glc
    x0[-1] = 0.1  # glx
    y = smallbone.compute_states(time_points=t, x0=x0)

    y = y[[0, 1, 2, 3, 4, 5]]
    for i in range(y.shape[0]):
        plt.plot(t, y[i, :])
    plt.legend(['glc', 'g1p', 'g6p', 'trh', 't6p', 'udg'])
    plt.title("glx=0.1")
    plt.show()
