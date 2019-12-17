import libsbml
import numpy as np
from math import exp, pow, log, log10, log2
# from bs4 import BeautifulSoup as bs
from pyADAPT import Model
from pyADAPT.bio import Reaction


class SBMLModel(Model):
    """The main purpose of this class is to support the yeast trehalose model from
    Smallbone et al. For other Model of different structure, some further extending
    might be needed.
    """
    def __init__(self, sbml_path, time_range=[]):
        self.sbml:libsbml.SBMLDocument = libsbml.readSBML(sbml_path)
        self.model:libsbml.Model = self.sbml.getModel()
        self.name = self.model.name
        self.notes = self.model.notes_string

        self.__reaction_list = list(self.model.reactions)
        self.reaction_list = list()
        self.species_list = list(self.model.species)
        # self.species_list = list()

        # stoichiometry matrix: row (effect, substrate), columns (cause, reaction)
        self.add_predictor(name='t', value=time_range)  # TODO move predictor to ADAPT
        
        self.stoich_matrix = np.zeros((len(self.model.species),
                                        len(self.model.reactions)))
        self.buildMatrix()
        self.buildReactions()
        self.buildStates()
        self.buildContext()
        self.initAssign()
        
        super().__init__()

    def buildMatrix(self):
        for r in self.model.reactions:
            # ignore those boundary species
            r_index = self.getReactionIndex(r)
            for sr in r.reactants:  # species reference of reactants
                s = self.model.getSpecies(sr.species)
                if not s.boundary_condition:
                    s_index = self.getSpeciesIndex(s)
                    self.stoich_matrix[s_index, r_index] = - sr.stoichiometry
            for sr in r.products:  # species reference of products
                s = self.model.getSpecies(sr.species)
                if not s.boundary_condition:
                    s_index = self.getSpeciesIndex(s)
                    self.stoich_matrix[s_index, r_index] = sr.stoichiometry

    def buildReactions(self):
        # convert libsbml.Reaction into pyADAPT.bio.Reaction, in a specific order
        for r in self.__reaction_list:
            # r_index = self.reaction_list.index(r)
            self.reaction_list.append(Reaction(r))
    
    def buildStates(self):
        # Following the convention of pyADAPT, the concentration of each species is a state
        for s in self.species_list:
            self.add_state(name=s.id, init=s.initial_concentration)

    def buildContext(self):
        context = {'log': log, 'log10': log10, 'log2': log2, 'exp': exp}
        for c in self.model.compartments:
            context[c.id] = c.size
        for p in self.model.parameters:
            context[p.id] = p.value
        for s in self.model.species:
            context[s.id] = s.initial_concentration
        # context.update(globals())
        context.update(self.states)
        self.context = context
        return context

    def getReactionIndex(self, r):
        return self.__reaction_list.index(r)
    
    def getSpeciesIndex(self, s):
        return self.species_list.index(s)
    
    
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
        # TODO they seems to be useless right now



    def odefunc(self, t, x, p):
        """
        Q: why odefunc take x (state) as an argument which is know to self
        A: the meaning of odefunc lies in setting the math form of the differential equation system, and the ability 
           to calculate the differentials at any time and states concentrations. If states is calculated from self, the
           ode function can only return one result.
        """
        v = self.reactions(t, x, p)
        dxdt = self.stoich_matrix.dot(v)
        return dxdt

    def reactions(self, t, x, p):
        # evaluate each reaction's flux(rate)
        self.updateContext(x)
        return [r.eval(self.context) for r in self.reaction_list]
        

    def updateContext(self, x):
        """ using solve_ivp, odefunc must take the states as a list rather than dictionary. The code should guarantee 
        thatn x is in the same order as species list
        """
        assert len(x) == len(self.species_list)
        for i in range(len(self.species_list)):
            s = self.species_list[i].id
            self.context[s] = x[i]


    def inputs(self, t):
        return super().inputs(t)

if __name__ == "__main__":
    from pprint import pprint
    smallbone = SBMLModel(r"data\trehalose\BIOMD0000000380_url.xml")
    print(smallbone.stoich_matrix)
    # pprint(smallbone.reaction_list)

    # print('context:', smallbone.context)
    x = np.ones(16)
    t_eval = np.linspace(0, 10)
    y = smallbone.compute_states([0, 10], x, t_eval=t_eval)
    import matplotlib.pyplot as plt
    for i in range(y.shape[0]):
        plt.plot(t_eval, y[i,:])
    names = [x.name for x in smallbone.species_list]
    plt.legend(names)
    plt.show()
