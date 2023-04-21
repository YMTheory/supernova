import numpy as np
from xsec.reaction_XS import reaction_XS

class NuP_XS(reaction_XS):

    def __init__(self) -> None:
        super().__init__(name="NuP")

    def diffXS(self, Ev, T, flavour):
        """
        input: T unit MeV, E unit MeV
        return: diffXS (cm2/MeV)
        """
        mp = 938 # MeV
        Emin = np.sqrt((1+T/2/mp) * (mp*T/2)) + T/2.
        if Ev < Emin:
            return 0
        diffxs = (1 + 466*T/Ev**2) * 4.83*1e-42
        return diffxs


    def totXS(self, Ev, flavour):
        """
        input: E unit MeV
        return: totXS (cm2/MeV)
        """
        mp = 938
        totxs = (2./mp+233*4/mp**2) * Ev**2 *4.83e-42
        return totxs

