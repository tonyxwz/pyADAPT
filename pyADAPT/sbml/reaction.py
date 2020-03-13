import math
from copy import deepcopy
import re

import libsbml
import numpy as np
from numexpr import evaluate
from pyADAPT.sbml.components import BaseNode


class Reaction(BaseNode):
    def __init__(self, reaction: libsbml.Reaction):
        super().__init__(reaction)
        # references, not actual species
        self.reactants = list(reaction.getListOfReactants())
        self.products = list(reaction.getListOfProducts())
        # inhibitors / activators
        self.modifiers = list(reaction.getListOfModifiers())
        self.reversible: bool = reaction.getReversible()

        self.kl = reaction.getKineticLaw()
        self.unit = libsbml.UnitDefinition.printUnits(
            self.kl.getDerivedUnitDefinition())
        self.formula = self.format_formula(self.kl.formula)

    def format_formula(self, text):
        # change pow(a, b) to a ** b
        regex1 = re.compile(r'pow\s*\(\s*\b(\w*)\b\s*,\s*\b(\w*)\b\s*\)')
        text = regex1.sub("(\\1 ** \\2)", text)
        # prepend reaction id_ to parameters id
        regex2 = re.compile(
            "(" + "|".join([p.id
                            for p in self.kl.getListOfParameters()]) + ")",
            re.VERBOSE,
        )
        text = regex2.sub(f"{self.id}_\\1", text)
        return text

    def compute_flux(self, context={}):
        """
        Calculate the fluxes of the reaction

        :param context: The required substrate concertations and parameters to
            calculate the flux.
        :return: flux (Float)
        """
        return evaluate(self.formula, context)

    def get_ce(self):
        # ce: chemical equation
        lhs = " + ".join(
            [str(x.stoichiometry) + " " + x.species for x in self.reactants])
        rhs = " + ".join(
            [str(x.stoichiometry) + " " + x.species for x in self.products])

        if self.reversible:
            arrow = " <=> "
        else:
            arrow = " -> "

        ce = arrow.join([lhs, rhs])
        if len(self.modifiers):
            ce += "; " + ",".join([x.species for x in self.modifiers])
        return ce

    def __repr__(self):
        s = f"Reaction {self.id}({self.name}): {self.get_ce()}"
        return s


if __name__ == "__main__":
    sbml = libsbml.readSBML("data/trehalose/smallbone.xml")
    model = sbml.getModel()
    context = dict()
    for c in model.getListOfCompartments():
        context[c.id] = c.size
    for p in model.getListOfParameters():
        context[p.id] = p.value
    for s in model.getListOfSpecies():
        context[s.id] = s.initial_concentration
    for r_ in model.getListOfReactions():
        r = Reaction(r_)
        print(r.name, ":", r, r.unit)
        print(r.id, ":", r.formula)
        for p in r.kl.getListOfParameters():
            print((p.id, p.value), end=" ")
        print("\n")
