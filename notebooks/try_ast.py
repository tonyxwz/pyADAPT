# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:51:02 2020

@author: 20175506
"""

from asteval import Interpreter
from math import log, log10

interpreter = Interpreter(usersyms={"glc": 100, "glc_0": 45})
glc_change = "log10(45/glc)"

# glc_change = interpreter("log10(45/glc)")
print(interpreter(glc_change))
