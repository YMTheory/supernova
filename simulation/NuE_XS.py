import numpy as np

class NuE_XS:

    def __init__(self) -> None:
        pass

    def diffXS(self, T, E, flavorID):
        me = 0.511 # MeV
        if E < (T/2+np.sqrt(T*(T+me))/2):
            return 0

        diffxs = 0
        sin2_thetaW = 0.23
        indexG = 1.732e-44

        if flavorID == 2 or flavorID == 4:
            epsi_p = -1 * sin2_thetaW
            epsi_m = 0.5 - sin2_thetaW
            diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))
            return diffxs

        if flavorID == 3 or flavorID == 5:
            epsi_p = -1*sin2_thetaW
            epsi_m = 0.5-sin2_thetaW
            diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))
            return diffxs
        
        if flavorID == 0:
            epsi_p = -1*sin2_thetaW
            epsi_m = -0.5-sin2_thetaW
            diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))
            return diffxs

        if flavorID == 1:
            epsi_p = -1*sin2_thetaW
            epsi_m = -0.5-sin2_thetaW
            diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))
            return diffxs

        return 0

    def redefine_ufunc(self):
        self.diffXS_ufunc = np.frompyfunc(self.diffXS, 3, 1)
