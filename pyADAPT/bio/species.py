class Species(object):
    def __init__(self, sp):
        self.id = sp.id
        self.name = sp.name
        self.conc = sp.initial_concentration
        self.compartment = sp.compartment