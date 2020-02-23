"""
The parameters that could be potentially changed during a simulation might
have two sources, the parameters in the model "model > listOfParameters"
and the parameters under a reaction "reaction > kineticLaw > 
listOfParameters".  

this can be problematic if the same name is defined in many reactions such
as "Vmax". In SBMl, the parameters have a scoop to operate. In pyADAPT,
the fluxes of a reaction is calculated with a context.

A proposed syntax for specifying the parent of the parameters
when calling from the command line is:

```sh
python -m pyADAPT analysis -f trehalose.xml -p R:hxt/Vmax
```

Where `hxt` is the name of the reaction and `Vmax` is the parameter that
the use want to "ADAPT"
"""

import libsbml
import numpy as np
from math import exp, pow, log, log10, log2
from cached_property import cached_property
from collections import OrderedDict

from pyADAPT.basemodel import BaseModel
from pyADAPT.bio.wrappers import Species, Compartment
from pyADAPT.bio.reaction import Reaction


class SBMLModel(BaseModel):
    """The main purpose of this class is to support the yeast trehalose model from
    Smallbone et al. For other Model of different structure, some further extending
    might be needed.
    """
    def __init__(self, sbml_path, time_range=[]):
        self.sbml: libsbml.SBMLDocument = libsbml.readSBML(sbml_path)
        self.model: libsbml.Model = self.sbml.getModel()
        assert self.model is not None
        self.name = self.model.name
        self.notes = self.model.notes_string
        # TODO, think I don't need the `sbml_reaction_list`, just convert em
        self.sbml_reaction_list = list(self.model.reactions)
        # TODO also here
        self.sbml_species_list = list(self.model.species)

        self.context = {"log": log, "log10": log10, "log2": log2, "exp": exp}
        # ! DO it
        self.reactions = OrderedDict()
        self.species = OrderedDict()
        # compartment sizes are constant for sure
        self.compartments = []
        for c_ in self.model.compartments:
            c = Compartment(c_)
            self.compartments.append(c)
            self.context[c.id] = c

        # self.parameters = list()
        for p_ in self.model.parameters:
            self.context[p_.id] = p_.value

        for s in self.model.species:
            self.context[s.id] = s.initial_concentration

        for s in self.sbml_species_list:
            self.add_state(name=s.id, init=s.initial_concentration)
        # self.buildStates()
        # self.context.update(self.states)
        self.initAssign()

        self.add_predictor(name="t", value=time_range)
        super().__init__()

    @cached_property
    def stoich_matrix(self):
        """
        stoichiometry matrix: row (effect or substrate)
                              columns (cause or reaction)
        """
        stoich_matrix = np.zeros(
            (len(self.model.species), len(self.model.reactions)))
        for r in self.model.reactions:
            # ignore those boundary species
            r_index = self.getReactionIndex(r)
            for sr in r.reactants:  # species reference of reactants
                s = self.model.getSpecies(sr.species)
                if not s.boundary_condition:
                    s_index = self.getSpeciesIndex(s)
                    stoich_matrix[s_index, r_index] = -sr.stoichiometry
            for sr in r.products:  # species reference of products
                s = self.model.getSpecies(sr.species)
                if not s.boundary_condition:
                    s_index = self.getSpeciesIndex(s)
                    stoich_matrix[s_index, r_index] = sr.stoichiometry
        return stoich_matrix

    @cached_property
    def reaction_list(self):
        # convert libsbml.Reaction into pyADAPT.bio.Reaction, in a specific order
        # TODO don't use this fancy cached property, I need to update the paramters in the reactions during the simulation, which means I need them all to be pyADAPT.bio.reaction.Reaction
        rl = list()
        for r in self.sbml_reaction_list:
            # r_index = self.reaction_list.index(r)
            rl.append(Reaction(r))
        return rl

    @cached_property
    def species_list(self):
        sl = list()
        for s in self.sbml_species_list:
            sl.append(Species(s))
        return sl

    def buildStates(self):
        # Following the convention of pyADAPT, the concentration of each species is a state
        # TODO observables in ADAPT convention is not taken into consideration
        #      maybe check whether a state exists in the dataset or not, then set
        #      observable flags accordingly
        for s in self.sbml_species_list:
            self.add_state(name=s.id, init=s.initial_concentration)

    def getReactionIndex(self, r):
        return self.sbml_reaction_list.index(r)

    def getSpeciesIndex(self, s):
        return self.sbml_species_list.index(s)

    def initAssign(self):
        """libSBML.InitialAssignment
        SBMl model support initial assignment to a symbol (species initial conc)
        This method is called after assign the value attribute of species, making it
        the higher priority in initialization.

        A initial assignment assigns the value of a formula comprised of global parameters
        to an entity such as substrate (concentration).
        A initial assignment can assign to a constant (default property of a parameter)
        """
        for ia in self.model.getListOfInitialAssignments():
            symbol = ia.symbol
            formula = libsbml.formulaToString(ia.math)
            self.states[symbol] = eval(formula, {}, self.context)

    def rulesToLambda(self):
        """libSBML.AssignmentRule
        the definition of Rule in SBML is identical to the definition of reactions
        in pyADAPT, libSBML::rule -> PyADAPT::Reactions
        """
        # TODO rules seems to be useless right now
        pass

    def get_unit(self, x):
        """ forget about the units for now
        they are only useful when i need to plot
        """
        pass

    def odefunc(self, t, x, p):
        """
        Q: why odefunc take x (state) as an argument which is know to self
        A: the meaning of odefunc lies in setting the math form of the differential
        equation system, and the ability to calculate the differentials at any time
        and states concentrations. If states is calculated from self, the ode
        function can only return one result.
        """
        v = self.fluxes(t, x, p)
        dxdt = self.stoich_matrix.dot(v)
        return dxdt

    def fluxes(self, t, x, p):
        # evaluate each reaction's flux(rate)
        self.updateContext(x)
        v = list()
        for r in self.reaction_list:
            # TODO if you want the parameters to be ADAPT, do it here
            v.append(r.compute_flux(self.context))
        return np.array(v)

    def updateContext(self, x):
        """ To use solve_ivp, odefunc must take the states as a list rather than
        dictionary. The code should guarantee that x is in the same order as species
        list.
        """
        # TODO this should be "update species"
        assert len(x) == len(self.sbml_species_list)
        for i in range(len(self.sbml_species_list)):
            s = self.sbml_species_list[i].id
            self.context[s] = x[i]

    def inputs(self, t):
        return super().inputs(t)


if __name__ == "__main__":
    from pprint import pprint

    smallbone = SBMLModel("data/trehalose/BIOMD0000000380_url.xml")
    print(smallbone.stoich_matrix)
    x = list()
    for s in smallbone.model.getListOfSpecies():
        x.append(s.initial_concentration)
    x = np.array(x)
    x[15] = 0.1

    # pprint(smallbone.reaction_list)

    # print('context:', smallbone.context)
    # x = np.random.rand(16) * 2
    # x = np.ones(16)
    t_eval = np.linspace(0, 10, 1000)
    y = smallbone.compute_states([0, 10], x, t_eval=t_eval)

    import matplotlib.pyplot as plt

    for i in range(y.shape[0]):
        plt.plot(t_eval, y[i, :], "--")
    names = [x.name for x in smallbone.sbml_species_list]
    plt.legend(names)
    plt.show()
