"""
Module sbml_model
===================

This module contains a convert to import a SBML model into the ADAPT model.

How to process SBML components:

- Model:
    - notes
    - annotation
    - list of compartments -> SBMlModel.parameters
        append to list of constants, the value is the size
    - list of inital assignments
        use in `__init__` to assign values to the species concentrations (states)
    - list of parameters -> SBMLModel.parameters
        these parameters are not the parameters in ADAPT. For example parameter
        `heat` is a constant in the smallbone model, while parameter `glc_0` is
        also a constant, but it is used in `initialassignment` to initialize
        species concentration. (for my project, the trehalose model)
    - list of reactions -> reaction.flux
        reactions are different from the ones in ADAPT because they contains
        local scope parameters, for example, `Vmax` exists in almost every
        reaction in the trehalose model. To import SBML model into pyADAPT,
        the parameters need to be renamed so that they don't conflict in the
        global scope. so, Vmax from hxt will become hxt_Vmax in the
        `SBMLModel.parameters`. Another to also store parameter under different
        scopes in this SBMLModel class.
    - list of rules -> rule.flux
        rules are similar to a property of the model, it should return a value
        calculated using the values known to the model already. They will be
        implemented as a reaction in `pyADAPT.SBMLModel`, but with an attribute
        marking they are actually rules not reactions.
    - list of species -> states
        Species are just the states in ADAPT models.
    - list of unit definitions
        fancy function for new students to work on.

The parameters that could be potentially changed during a simulation have two
sources, the parameters in the model "model > listOfParameters"
and the parameters in a reaction "reaction > kineticLaw > listOfParameters".

- "model > listOfParameters" are just some initial conditions
- "reaction > kineticLaw > listOfParameters" are the parameters that truely matters

A proposed syntax for specifying the parent of the parameters
when calling from the command line is:

```sh
python -m pyADAPT analysis -f trehalose.xml -p R:hxt/Vmax
```

Where `hxt` is the name of the reaction and `Vmax` is the parameter that
the user want to "ADAPT"
"""

from collections import OrderedDict
from math import exp, log, log2, log10, pow

import libsbml
import numpy as np
import tempfile  # TODO return Model object directly

from pyADAPT.basemodel import BaseModel
from pyADAPT.sbml.reaction import Reaction
from pyADAPT.sbml.components import Compartment, Species, AssignmentRule, SParameter


class SBMLModel(BaseModel):
    """ The main purpose of this class is to support the yeast trehalose model from
    Smallbone et al. For other Model of different structure, further extending
    might be required.
    """

    def __init__(self, sbml_path):
        """Initialize a ADAPT model from a sbml file

        Parameters
        ----------
        sbml_path : String
            The path to the sbml model 'xml' file
        """
        self.sbml: libsbml.SBMLDocument = libsbml.readSBML(sbml_path)
        # self.sbml_model: libsbml.Model = self.sbml.getModel()
        assert self.sbml.model is not None
        self.name = self.sbml.model.name
        self.notes = self.sbml.model.notes_string

        self.math_functions = {"log": log, "log10": log10, "log2": log2, "exp": exp}
        # add to _parameters
        for c in self.sbml.getModel().getListOfCompartments():
            self.add_parameter(c.id, c.size, False, lb=0, ub=np.inf)

        # self.getRules(context)
        rules = self.sbml.model.getListOfRules()
        self.rules = OrderedDict()
        for r in rules:
            self.rules[r.variable] = AssignmentRule(r)

        for p in self.sbml.getModel().getListOfParameters():
            r = rules.get(p.id)
            if r is None:
                self.add_parameter(p.id, p.value, not p.constant, lb=0, ub=np.inf)

        # add reactions as pyADAPT.sbml.reaction.Reaction
        self.reactions = OrderedDict()
        for sbml_rxn in self.sbml.model.getListOfReactions():
            rxn = Reaction(sbml_rxn)
            self.reactions[rxn.id] = rxn
            # move parameters from reaction to model scope, by add prefix
            for param in rxn.kl.getListOfParameters():
                param_name = "_".join([rxn.id, param.id])
                self.add_parameter(
                    name=param_name,
                    value=param.value,
                    vary=False,
                    lb=0,
                    ub=np.inf
                    # parent=rxn.id,
                )

        # ias: list of initial assignments
        state_order = []
        initial_states = []
        flux_order = self.reactions.keys()

        ias = self.sbml.getModel().getListOfInitialAssignments()
        for s in self.sbml.getModel().getListOfSpecies():
            name = s.id
            ia = ias.get(name)
            if ia is not None:
                formula = libsbml.formulaToString(ia.getMath())
                conc = eval(formula, {}, {p[0]: p[1] for p in self._parameters})
            elif s.boundary_condition and s.initial_concentration == 0:
                conc = 1.0
            else:
                conc = s.initial_concentration

            if s.boundary_condition:
                self.add_parameter(s.id, conc, False)
            else:
                state_order.append(s.id)
                initial_states.append(conc)
            # self.add_state(name=s.id, value=conc, observable=True)
        super().__init__(
            state_order=state_order,
            initial_states=initial_states,
            flux_order=list(flux_order),
        )

        self.stoich_matrix = self.get_stoich_matrix()
        self.symbol_table = {}
        self.symbol_table.update(self.math_functions)

    def get_stoich_matrix(self):
        """
        stoichiometry matrix: row (effect or substrate/species)
                                columns (cause or reaction)
        """
        mat = np.zeros((len(self.state_order), len(self.reactions)))
        # maybe it's a good idea to always refer to the sbml for ordering
        species_list = self.state_order
        for j, r in enumerate(self.reactions.values()):
            for species_ref in r.reactants:
                reactant = self.sbml.model.species.get(species_ref.species)
                if not reactant.boundary_condition:
                    i = species_list.index(species_ref.species)
                    mat[i, j] = -species_ref.stoichiometry
            for species_ref in r.products:
                product = self.sbml.model.species.get(species_ref.species)
                if not product.boundary_condition:
                    i = species_list.index(species_ref.species)
                    mat[i, j] = species_ref.stoichiometry
        return mat

    @property
    def symbols(self):
        """parameters, states,"""
        table = self.math_functions
        table.update(self.parameters["value"].to_dict())
        table.update(self.states["value"].to_dict())
        return table

    def fluxes(self, t, x, p):
        # evaluate each reaction's flux(rate)
        # Note: it's a bit weird since all the fluxes are `eval`ed, `p` is not
        # used in this method.
        v = list()
        for r in self.reactions.values():
            v.append(r.compute_flux(self.symbols))
        v = np.array(v)
        return v

    def state_ode(self, t, x, p):
        """
        calculate the derivatives in biologist's style
        p: np.ndarray / pandas.Series
        """
        # assign p (parameters) and x (states)
        # self.states.loc[:, "value"] = x
        v = self.fluxes(t, x, p)
        dxdt = np.dot(self.stoich_matrix, v)
        return dxdt


def SBMLBuilder(sbml) -> SBMLModel:
    # TODO
    # builder pattern to avoid the adding libsbml object, which cannot be pickled
    # when using multiprocessing
    pass


if __name__ == "__main__":
    from pprint import pprint
    from copy import deepcopy

    smallbone = SBMLModel("data/trehalose/smallbone.xml")
    print(smallbone.state_order)
    print(smallbone.flux_order)
    print(smallbone.stoich_matrix)
    # x1 = deepcopy(smallbone.states["value"].values)
    # x0 = np.array(
    #     [
    #         0.09765,
    #         0.10000,
    #         2.67500,
    #         0.05000,
    #         0.02000,
    #         0.70000,
    #         1.28200,
    #         2.52500,
    #         1.00000,
    #         0.62500,
    #         1.00000,
    #         1.00000,
    #         0.28150,
    #         0.64910,
    #         1.00000,
    #         100.00000,
    #     ]
    # )
    # x0[15] = 0.1
    # x0[0] = 100
    # t_eval = np.linspace(0, 10, 1000)
    # # import cProfile
    # # cProfile.run("""smallbone.compute_states(new_params=0.5,
    # #                              time_points=t_eval,
    # #                              x0=x,
    # #                              new_param_names=['hxt_Vmax'])""")
    # y = smallbone.compute_states(time_points=t_eval, x0=x0)

    # # print(y)

    # import matplotlib.pyplot as plt

    # for i in range(y.shape[0]):
    #     plt.plot(t_eval, y[i, :], "--")
    # # names = [x.name for x in smallbone.states.values()]
    # print(smallbone.states["value"])
    # plt.legend(smallbone.states.name)
    # plt.show()
    for i in smallbone.parameters.index:
        print(i)

    for k, v in smallbone.reactions.items():
        print(k, "=", v.text_formula)

