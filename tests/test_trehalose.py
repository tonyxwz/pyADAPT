"""Test case goal
Test the smallbone2011 model
1. runnable
2. give same result as the matlab code by David
"""

import unittest
from pyADAPT.sbml.sbml_model import SBMLModel
import os
import sys


class TestSmallbone2011(unittest.TestCase):
    def setUp(self):
        sbml_path = os.path.join(os.path.dirname(__file__), "..", "data",
                                 "trehalose", "smallbone.xml")
        self.model = SBMLModel(sbml_path)

    def test_basic(self):
        self.assertIsInstance(self.model, SBMLModel)

    def test_constants(self):
        pass
