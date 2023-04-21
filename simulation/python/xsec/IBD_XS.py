import numpy as np

from xsec.reaction_XS import reaction_XS

class IBD_XS(reaction_XS) :
    def __init__(self) -> None:
        super().__init__(name="IBD")

    def diffXS(self, Ev, T, flavour):
        return 0.


    def totXS(self, Ev, flavour):
        """
        return IBD xs: unit cm2
        """
        Ethr = 1.8022
        if Ev < Ethr:
            return 0
        else:
            mn = 939.57
            mp = 938.27
            delta_mnp = mn - mp
            me = 0.511
            Ee = Ev - delta_mnp
            Pe = np.sqrt(Ee**2 - me**2)
            totxs = 0.0952 * (Ee*Pe) * 1e-42 
            return totxs
