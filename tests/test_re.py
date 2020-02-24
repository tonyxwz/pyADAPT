from pyADAPT.bio.sbml_model import SBMLModel
import unittest
import os

test_formulas = {
    "pgi": "cell * pow(pgi_shock, heat) * pgi_Vmax / pgi_Kg6p * (g6p - f6p / pgi_Keq) / (1 + g6p / pgi_Kg6p + f6p / pgi_Kf6p)",
    "hxt": "cell * pow(hxt_shock, heat) * hxt_Vmax * (glx - glc) / hxt_Kglc / (1 + (glx + glc) / hxt_Kglc + hxt_Ki * glx * glc / pow(hxt_Kglc, 2))",
    "hxk": "cell * pow(hxk_shock, heat) * hxk_Vmax / (hxk_Kglc * hxk_Katp) * (glc * atp - g6p * adp / hxk_Keq) / ((1 + glc / hxk_Kglc + g6p / hxk_Kg6p + t6p / hxk_Kit6p) * (1 + atp / hxk_Katp + adp / hxk_Kadp))",
    "pgm": "cell * pow(pgm_shock, heat) * pgm_Vmax / pgm_Kg6p * (g6p - g1p / pgm_Keq) / (1 + g6p / pgm_Kg6p + g1p / pgm_Kg1p)",
    "tpp": "cell * pow(tpp_shock, heat) * tpp_Vmax * t6p / tpp_Kt6p / (1 + t6p / tpp_Kt6p)",
    "tps": "cell * tps_activity * pow(tps_shock, heat) * tps_Vmax * g6p * udg / (tps_Kg6p * tps_Kudg) / ((1 + g6p / tps_Kg6p) * (1 + udg / tps_Kudg))",
    "nth": "cell * pow(nth_shock, heat) * nth_Vmax * trh / nth_Ktrh / (1 + trh / nth_Ktrh)",
    "ugp": "cell * pow(ugp_shock, heat) * ugp_Vmax * utp * g1p / (ugp_Kutp * ugp_Kg1p) / (ugp_Kiutp / ugp_Kutp + utp / ugp_Kutp + g1p / ugp_Kg1p + utp * g1p / (ugp_Kutp * ugp_Kg1p) + ugp_Kiutp / ugp_Kutp * udg / ugp_Kiudg + g1p * udg / (ugp_Kg1p * ugp_Kiudg))",
}


cwd = os.path.dirname(__file__)
xml = os.path.join(cwd, "..", "data/trehalose/smallbone.xml")


class TestRegex(unittest.TestCase):
    def test_formula(self):
        # output = {}
        smallbone = SBMLModel(xml)
        for k, r in smallbone.reactions.items():
            self.assertEqual(test_formulas[k], r.uncompiled_formula)

