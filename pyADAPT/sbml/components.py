"""the SBML components we need to consider in this project
create a shadow class for every node

    - Reaction: compute_flux(context)
    - Compartment: cell=1, medium=1
    - Species:
    - UnitDefinition
    - Parameter
    - InitialAssignment
"""

import libsbml
from numexpr import evaluate


class BaseNode(object):
    def __init__(self, node: libsbml.SBase):
        self.id = node.id
        self.meta_id = node.meta_id
        self.name = node.name
        self.value: float

    def __add__(self, other):
        return self.value + other

    def __radd__(self, other):
        return other + self.value

    def __sub__(self, other):
        return self.value - other

    def __rsub__(self, other):
        return other - self.value

    def __mul__(self, other):
        return self.value * other

    def __rmul__(self, other):
        return other * self.value

    def __truediv__(self, other):
        return self.value / other

    def __rtruediv__(self, other):
        return other / self.value

    def __floordiv__(self, other):
        return self.value // other

    def __rfloordiv__(self, other):
        return other // self.value

    def __mod__(self, other):
        return self.value % other

    def __rmod__(self, other):
        return other % self.value

    def __pow__(self, other):
        return self.value**other

    def __rpow__(self, other):
        return other**self.value

    def __neg__(self):
        return -self.value

    def __pos__(self):
        return self.value

    def __abs__(self):
        return abs(self.value)

    def __iadd__(self, other):
        self.value += other

    def __isub__(self, other):
        self.value -= other

    def __eq__(self, other):
        return self.value == other

    def __ne__(self, other):
        return self.value != other

    def __lt__(self, other):
        return self.value < other

    def __le__(self, other):
        return self.value <= other

    def __gt__(self, other):
        return self.value > other

    def __ge__(self, other):
        return self.value >= other

    def __repr__(self):
        return str((self.name, self.id, self.value))


class Compartment(BaseNode):
    def __init__(self, comp: libsbml.Compartment):
        super(Compartment, self).__init__(comp)
        self.outside = comp.getOutside()
        self.value = comp.size

    @property
    def size(self):
        return self.value


class Species(BaseNode):
    def __init__(self, sp: libsbml.Species):
        super().__init__(sp)
        self.value = sp.initial_concentration
        self.compartment = sp.compartment
        self.boundary_condition = sp.getBoundaryCondition()

    @property
    def initial_concentration(self):
        return self.value


class AssignmentRule(BaseNode):
    def __init__(self, rule: libsbml.AssignmentRule):
        super().__init__(rule)
        self.variable = rule.getVariable()  # target variable
        self.formula = libsbml.formulaToString(rule.math)

    def __repr__(self):
        return f"pyADAPT AssginmentRule: {self.variable} = {self.formula}"

    def get_value(self, context):
        return evaluate(self.formula, {}, context)


class SParameter(BaseNode):
    def __init__(self, pnode: libsbml.Parameter):
        """ 'S' for SBML
        http://sbml.org/Software/libSBML/docs/python-api/classlibsbml_1_1_parameter.html
        By default, a SBML parameter
        A wrapper over sbml parameters to provide the support for assignment rule
        """
        super().__init__(pnode)
        self.units = pnode.units


if __name__ == "__main__":
    sbml = libsbml.readSBML("data/trehalose/smallbone.xml")
    model = sbml.model
    compartments = model.getListOfCompartments()
    # print(compartments[0], compartments[1])
    cell = Compartment(compartments.getElementBySId("cell"))

    print(
        cell + 0.5 == cell.value + 0.5,
        0.5 + cell == 0.5 + cell.value,
        cell - 0.5 == cell.value - 0.5,
        0.5 - cell == 0.5 - cell.value,
        cell * 4.5 == cell * 4.5,
        4.5 * cell == 4.5 * cell.value,
        cell / 8.3 == cell.value / 8.3,
        8.3 / cell == 8.3 / cell.value,
        cell // 8.3 == cell.value // 8.3,
        8.3 // cell == 8.3 // cell.value,
        cell % 4 == cell.value % 4,
        4 % cell == 4 % cell.value,
        cell**4 == cell.value**4,
        4**cell == 4**cell.value,
    )

    medium = Compartment(compartments.getElementBySId("medium"))
    print(
        cell == medium,
        cell > medium,
        cell < medium,
        cell <= medium,
        cell >= medium,
        cell != medium,
    )
    print(cell, medium)
    for ar in model.getListOfRules():
        print(AssignmentRule(ar))
