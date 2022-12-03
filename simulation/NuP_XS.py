import numpy as np

class NuP_XS:

    def __init__(self) -> None:
        pass

    def diffXS(self, T, E):
        """
        input: T unit MeV, E unit MeV
        return: diffXS (cm2/MeV)
        """
        mp = 938 # MeV
        Emin = np.sqrt((1+T/2/mp) * (mp*T/2)) + T/2.
        if E < Emin:
            return 0
        diffxs = (1 + 466*T/E**2) * 4.83*1e-42
        return diffxs


    def totXS(E):
        """
        input: E unit MeV
        return: totXS (cm2/MeV)
        """
        mp = 938
        totxs = (2./mp+233*4/mp**2) * E**2 *4.83e-42
        return totxs
