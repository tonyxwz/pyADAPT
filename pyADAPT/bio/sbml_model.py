"""
Module `sbml_model`
===================

This module offers an ADAPT model that can be converted from a SBML model.

The parameters that could be potentially changed during a simulation have two
sources, the parameters in the model "model > listOfParameters"
and the parameters in a reaction "reaction > kineticLaw > listOfParameters".

- "model > listOfParameters" are just some initial conditions
- "reaction > kineticLaw > listOfParameters" are the parameters that truely matters

this can be problematic if the same name , such as "Vmax", is defined in many
reactions. In SBMl, the parameters have a scope to operate. However, ADAPT
controls the parameters only from a global scope. A function to pick out the
parameters is required.

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
    def __init__(self, sbml_path, time_range=[], adapt_params=[]):
        """Initialize a ADAPT model from a sbml file
        
        Parameters
        ----------
        sbml_path : String
            The path to the sbml model 'xml' file
        time_range : array(2), optional
            array([t0, tf]), by default []
        adapt_params : list of strings, optional
            the parameter that need to be "adapted", for example,
            ['glc_0', 'hxt/Vmax'], by default []
        """
        self.sbml: libsbml.SBMLDocument = libsbml.readSBML(sbml_path)
        self.sbml_model: libsbml.Model = self.sbml.getModel()
        assert self.sbml_model is not None

        # Start to parse the sbml into an ADAPT model
        self.name = self.sbml_model.name
        self.notes = self.sbml_model.notes_string

        # math functions in the sbml formulas
        self.math_functions = {
            "log": log,
            "log10": log10,
            "log2": log2,
            "exp": exp
        }
        # TODO add parameters (dummy)
        # rules -> sbml.parameter
        rules = self.sbml_model.getListOfRules()
        for p in self.sbml_model.getListOfParameters():
            rule_tmp = rules.get(p.id)
            if rule_tmp is None:
                # constant parameter
                self.constants[p.id] = p.value

        # ? consider use BaseModel.constants
        # self.compartments = OrderedDict()
        for c in self.sbml_model.compartments:
            self.constants[c.id] = c.size
            # self.compartments[c.id] = Compartment(c)

        # add reactions which are ordered and whose states are updatable
        self.reactions = OrderedDict()
        for r in self.sbml_model.getListOfReactions():
            self.reactions[r.id] = Reaction(r)

        # initial_assignment -> species.initial_concentration
        self.species = OrderedDict()  # remember the order
        ia_list = self.sbml_model.getListOfInitialAssignments()
        for s in self.sbml_model.getListOfSpecies():
            sp = Species(s)
            ia = ia_list.get(sp.id)
            if ia is not None:
                formula = compile(libsbml.formulaToString(ia.math), "<string>",
                                  "eval")
                init_conc = eval(formula, {}, self.context)
            else:
                init_conc = s.initial_concentration
            self.species[s.id] = s
            self.add_state(name=s.id, init=init_conc)

        self.add_predictor(name="t", value=time_range)
        super().__init__()

    def add_species(self, s: libsbml.Species):
        pass

    def add_reaction(self, r: libsbml.Reaction):
        pass

    @cached_property
    def stoich_matrix(self):
        return self.get_stoich_matrix()

    def get_stoich_matrix(self):
        """
        stoichiometry matrix: row (effect or substrate/species)
                              columns (cause or reaction)
        """
        mat = np.zeros((len(self.species), len(self.reactions)))
        species_list = list(self.species)  # list of keys
        for j, r in enumerate(self.reactions.values()):
            for species_ref in r.reactants:
                reactant = self.species[species_ref.species]
                if not reactant.boundary_condition:
                    i = species_list.index(species_ref.species)
                    mat[i, j] = -species_ref.stoichiometry
            for species_ref in r.products:
                product = self.species[species_ref.species]
                if not product.boundary_condition:
                    i = species_list.index(species_ref.species)
                    mat[i, j] = species_ref.stoichiometry
        return mat

    @property
    def context(self):
        return self.get_context()

    def get_context(self):
        """
        TODO automatically update all the context
        - compartment sizes
        - species states concentrations
        - parameters (glc_0)
        - reaction parameters
        """
        context = {}
        context.update(self.constants)
        context.update(self.math_functions)
        context.update(self.parameters)
        context.update(self.states)
        return context

    def rulesToLambda(self):
        """libSBML.AssignmentRule
        the definition of Rule in SBML is identical to the definition of reactions
        in pyADAPT, libSBML::rule -> PyADAPT::Reactions
        glc_change = log10(glc / glc_0) at any time
        self.rules['glc_change'] = compile('log10(glc/glc_0)', '<string>', 'eval')
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
        # print(x)
        self.update_states(x)
        v = list()
        for r in self.reactions.values():
            # TODO if you want the parameters to be ADAPT, do it here
            v.append(r.compute_flux(self.context))
        return np.array(v)

    def pick_params_for_reaction(self, p):
        pass

    def update_states(self, x):
        """ To use solve_ivp, odefunc must take the states as a list rather than
        dictionary. The code should guarantee that x is in the same order as species
        list.
        """
        # TODO this should be "update species"
        assert len(x) == len(self.species)
        for species_id, new_state in zip(self.species.keys(), x):
            self.states[species_id] = new_state

    def inputs(self, t):
        # Not supported
        return super().inputs(t)


if __name__ == "__main__":
    from pprint import pprint

    smallbone = SBMLModel("data/trehalose/BIOMD0000000380_url.xml")
    print(smallbone.stoich_matrix)
    x = list()
    for s in smallbone.species.values():
        x.append(s.initial_concentration)
    x = np.array(x)
    x[15] = 0.1

    t_eval = np.linspace(0, 10, 100)
    y = smallbone.compute_states([0, 10], x, t_eval=t_eval)

    import matplotlib.pyplot as plt

    for i in range(y.shape[0]):
        plt.plot(t_eval, y[i, :], "--")
    names = [x.name for x in smallbone.species.values()]
    plt.legend(names)
    plt.show()
