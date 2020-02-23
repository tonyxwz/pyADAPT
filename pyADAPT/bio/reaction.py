import math
from copy import deepcopy

import libsbml
import numpy as np

from pyADAPT.bio.wrappers import BaseNode, Compartment, Species


class Reaction(BaseNode):
    def __init__(self, reaction: libsbml.Reaction):
        super().__init__(reaction)
        self.reactants = list(reaction.getListOfReactants())
        self.products = list(reaction.getListOfProducts())
        # inhibitors / activators
        self.modifiers = list(reaction.getListOfModifiers())
        self.reversible: bool = reaction.getReversible()

        self.kl = reaction.getKineticLaw()
        self.unit = libsbml.UnitDefinition.printUnits(
            self.kl.getDerivedUnitDefinition())
        self.formula = compile(self.kl.formula, "<string>", "eval")

        # TODO context need to be updated during the simulation
        self.context = dict()
        for p in self.kl.getListOfParameters():
            self.context[p.id] = p.value
        # TODO: listofspecies in Reaction actually stores the reference

    def compute_flux(self, context={}):
        """
        Calculate the fluxes of the reaction
        
        :param context: The required substrate concertations and parameters to
            calculate the flux.
        :return: flux (Float)
        """
        # TODO find out is `deepcopy` necessary?
        # context = deepcopy(context)
        self.context.update(context)  # update the context
        return eval(self.formula, {}, self.context)

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
    sbml = libsbml.readSBML("data/trehalose/BIOMD0000000380_url.xml")
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
        print(r.name, ":", r, r.compute_flux(context))
