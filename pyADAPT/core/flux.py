import numpy as np
import re


class Flux(object):
    """Flux as denoted by vector `f` in the literatures
    the flow from one species to another
    """
    def __init__(self, xs, xe, equation):
        """ A flux aka reaction in the ADAPT paper and Tiemann's thesis, is a
        mathematical expression used in the model. Sometimes they are real
        measurable values. Sometimes they are just some helper values to ease
        definition of ODEs.
        
        Parameters
        ----------
        `xs` : start species
        `xe` : end species
        `value` : flux
        """
        self.xs = xs
        self.xe = xe
        self.equation = equation # p_1 * p_2
