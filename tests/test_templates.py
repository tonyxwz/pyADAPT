import unittest
from pyADAPT.convert import Converter


class Test_Templates(unittest.TestCase):
    def setUp(self):
        self.sbml_path = "data/trehalose/smallbone.xml"
        self.converter = Converter(self.sbml_path)

    def test_converted_model(self):
        pass


tmplt.render(
    {
        "model_name": "TestModel",
        "model_description": "Model description",
        "predictors": [{"name": "g2", "value": 63}, {"name": "t4", "value": 56}],
        "constants": [{"name": "u1", "value": 1}, {"name": "u2", "value": 3}],
        "parameters": [{"name": "k1", "value": 4, "vary": True, "lb": 0}],
        "states": [{"name": "s1", "init": 1}, {"name": "s2", "init": 2}],
    }
)

