"""TTM: Tony's Toy Model
This model is created with the intention to provide an easier tutorial model other than the default toy model created in
2013 by Natal van Riel et al.
"""

import numpy as np
import cobra
from pyADAPT import Model


class TonysToyModel(Model):
    """ 
    # Tony's Toy Model
    A kinetic model of Trehalose biosynthesis in Saccharomyces Cerevisiae
    
    # Trehalose cycle

    # References
    https://www.ncbi.nlm.nih.gov/pubmed/21943906
    """

    def __init__(self):
        """substrates:
        1. GLC: cellular glucose
        2. GLX: medium glucose
        3. G6P: glucose 6-phosphate
        4. F6P: fructose 6-phophate
        5. G1P: glucose 1-phosphate
        6. UTP: uridine triphosphate
        7. UDP: uridine diphosphate
        8. UDG: UDP-glucose
        9. PPi: diphosphate (P2O7^4-)
        10. Pi: phosphate (PO4^3-)
        11. T6P: Trehalose 6-phosphate
        12. TRH: trehalose
        13. H+: Hydrogen ion (proton)
        14. ATP: Adenosine triphosphate
        15. ADP: Adenosine diphosphate
        16. H2O: water

        Enzymes:

        """
        self.name = 'Tony\'s toy model'
        self.description = self.__doc__
        self.add_predictor(name="time", value=[0, 10])
        
        #
        super().__init__()


    def reactions(self, t, x, p):
        """
        reactions:
        1. Hexokinase: GLC + ATP -> G6P + ADP + H+
        2. Phosphoglucomutase: G6P <=> G1P
        3. UDP-glucose phosphorylase: G1P + UTP + H+ -> UDG + PPi
        4. T6P synthase: G6P + UDG -> T6P + UDP + H+
        5. T6P + H2O -> TRH + Pi
        6. TRH + H2O -> 2 GLC
        7. GLX <=> GLC
        8. G6P <=> F6P
        """
        GLC = x['GLC']
        GLX = x['GLX']
        G6P = x['G6P']
        F6P = x['F6P']
        G1P = x['G1P']
        UTP = x['UTP']
        UDP = x['UDP']
        UDG = x['UDG']
        PPi = x['PPi']
        Pi  = x['Pi']
        T6P = x['T6P']
        TRH = x['TRH']
        H1  = x['H1']
        ATP = x['ATP']
        ADP = x['ADP']
        H2O = x['H2O']

        Vcell = 1  # volume of the cell

        

        
        v_hexokinase = V_max_hexokinase * d

    def input(self, t):
        """ TODO
        """