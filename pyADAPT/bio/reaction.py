import libsbml
import numpy as np

from copy import deepcopy
import math

class Reaction(object):
    def __init__(self, rnode:libsbml.Reaction):
        self.id = rnode.getId()
        self.name = rnode.getName()
        self.reactants = rnode.getListOfReactants()
        self.products = rnode.getListOfProducts()
        self.modifiers = rnode.getListOfModifiers()  # inhibitors / activators
        self.reversible = rnode.getReversible()

        self.kl = rnode.getKineticLaw()
        self.formula = compile(self.kl.formula, '<string>', 'eval')

        self.context = dict()
        for p in self.kl.getListOfParameters():
            self.context[p.id] = p.value

    def eval(self, context={}):
        # evaluate the flux of Reaction given all the concentrations
        context = deepcopy(context)
        context.update(self.context)
        return eval(self.formula, globals(), context)

    def get_ce(self):
        # ce: chemical equation
        lhs = " + ".join([str(x.stoichiometry)+" "+x.species for x in self.reactants])
        rhs = " + ".join([str(x.stoichiometry)+" "+x.species for x in self.products])

        if self.reversible:
            arrow = ' <=> '
        else:
            arrow = ' -> '

        ce = arrow.join([lhs, rhs])
        if len(self.modifiers):
            ce += '; ' + ','.join([x.species for x in self.modifiers])
        return ce

    def __repr__(self):
        # s = "{:>20}: {:>10}".format(f"Reaction {self.id}({self.name})", self.get_ce())
        s = f"Reaction {self.id}({self.name}): {self.get_ce()}"
        s = "{:<40} {}".format(f'Reaction {self.id}({self.name}):', self.get_ce())
        return s


if __name__ == "__main__":
    sbml = libsbml.readSBML(r"data\trehalose\BIOMD0000000380_url.xml")
    model = sbml.getModel()
    hxk = model.getReaction('pgm')
    hxk = Reaction(hxk)

    context = dict()
    for c in model.getListOfCompartments():
        context[c.id] = c.size
    for p in model.getListOfParameters():
        context[p.id] = p.value
    for s in model.getListOfSpecies():
        context[s.id] = s.initial_concentration
    print(hxk.eval(context))
    print(hxk)
