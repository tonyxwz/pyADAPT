import libsbml
import numpy as np
from math import exp, pow, log, log10, log2
from cached_property import cached_property

# from bs4 import BeautifulSoup as bs
from pyADAPT import Model
from pyADAPT.bio import Reaction, Species


class SBMLModel(Model):
    """The main purpose of this class is to support the yeast trehalose model from
    Smallbone et al. For other Model of different structure, some further extending
    might be needed.
    """

    def __init__(self, sbml_path, time_range=[]):
        self.sbml: libsbml.SBMLDocument = libsbml.readSBML(sbml_path)
        self.model: libsbml.Model = self.sbml.getModel()
        self.name = self.model.name
        self.notes = self.model.notes_string

        self.sbml_reaction_list = list(self.model.reactions)

        self.sbml_species_list = list(self.model.species)
        # self.species_list = list()

        # stoichiometry matrix: row (effect, substrate), columns (cause, reaction)
        # TODO move predictor to ADAPT
        self.add_predictor(name="t", value=time_range)
        self.context = {"log": log, "log10": log10, "log2": log2, "exp": exp}
        for c in self.model.compartments:
            self.context[c.id] = c.size
        for p in self.model.parameters:
            self.context[p.id] = p.value
        for s in self.model.species:
            self.context[s.id] = s.initial_concentration

        self.buildStates()
        # self.context.update(self.states)
        self.initAssign()

        super().__init__()

    @cached_property
    def stoich_matrix(self):
        stoich_matrix = np.zeros((len(self.model.species), len(self.model.reactions)))
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
            self.states[symbol] = eval(formula, globals(), self.context)

    def rulesToLambda(self):
        """libSBML.AssignmentRule
        the definition of Rule in SBML is identical to the definition of reactions
        in pyADAPT, libSBML::rule -> PyADAPT::Reactions
        """
        # TODO rules seems to be useless right now

    def odefunc(self, t, x, p):
        """
        Q: why odefunc take x (state) as an argument which is know to self
        A: the meaning of odefunc lies in setting the math form of the differential
        equation system, and the ability to calculate the differentials at any time
        and states concentrations. If states is calculated from self, the ode
        function can only return one result.
        """
        v = self.reactions(t, x, p)
        dxdt = self.stoich_matrix.dot(v)
        return dxdt

    def reactions(self, t, x, p):
        # evaluate each reaction's flux(rate)
        self.updateContext(x)
        return [r.eval(self.context) for r in self.reaction_list]

    def updateContext(self, x):
        """ To use solve_ivp, odefunc must take the states as a list rather than
        dictionary. The code should guarantee that x is in the same order as species
        list.
        """
        assert len(x) == len(self.sbml_species_list)
        for i in range(len(self.sbml_species_list)):
            s = self.sbml_species_list[i].id
            self.context[s] = x[i]

    def inputs(self, t):
        return super().inputs(t)


if __name__ == "__main__":
    from pprint import pprint

    smallbone = SBMLModel(r"data\trehalose\BIOMD0000000380_url.xml")
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
    t_eval = np.linspace(0, 10, 50)
    y = smallbone.compute_states([0, 10], x, t_eval=t_eval)
    import matplotlib.pyplot as plt

    for i in range(y.shape[0]):
        plt.plot(t_eval, y[i, :], "--")
    names = [x.name for x in smallbone.sbml_species_list]
    plt.legend(names)
    plt.show()
