import numpy as np

from xsec.reaction_XS import reaction_XS

class NuE_XS(reaction_XS):

    def __init__(self) -> None:
        super().__init__(name="NuE")

    def diffXS(self, Ev, T, flavorID):
        me = 0.511 # MeV
        if Ev < (T/2+np.sqrt(T*(T+me))/2):
            return 0

        diffxs = 0
        sin2_thetaW = 0.23
        indexG = 1.732e-44

        if flavorID == 2 or flavorID == 4:
            epsi_p = -1 * sin2_thetaW
            epsi_m = 0.5 - sin2_thetaW
            diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*np.power((1-T/Ev),2)-epsi_p*epsi_m*me*T/(Ev*Ev))
            return diffxs

        if flavorID == 3 or flavorID == 5:
            epsi_p = -1*sin2_thetaW
            epsi_m = 0.5-sin2_thetaW
            diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*np.power((1-T/Ev),2)-epsi_p*epsi_m*me*T/(Ev*Ev))
            return diffxs
        
        if flavorID == 0:
            epsi_p = -1*sin2_thetaW
            epsi_m = -0.5-sin2_thetaW
            diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*np.power((1-T/Ev),2)-epsi_p*epsi_m*me*T/(Ev*Ev))
            return diffxs

        if flavorID == 1:
            epsi_p = -1*sin2_thetaW
            epsi_m = -0.5-sin2_thetaW
            diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*np.power((1-T/Ev),2)-epsi_p*epsi_m*me*T/(Ev*Ev))
            return diffxs

        return 0


    def totXS(self, Ev, flavour): 
        return 0.
