import numpy as np
import scipy.integrate as integrate

class CEvNS_XS :
    
    def __init__(self, Z, N, mN):
        self.Z = Z
        self.N = N
        self.mN = mN # unit: eV
   
    def totXS_intgpart(self, theta, Ev, A):
        q = np.sqrt(2*Ev**2*(1-np.cos(theta)))
        return (1 + np.cos(theta)) * np.sin(theta) # firstly, omit form factor

    def totXS(self, Ev):
        """
        return total cross section with neutrino energy.
        Ev unit: MeV
        XS unit: cm2 
        """
        GF = 1.1664e-5
        A = self.N + self.Z
        eVminus1tocm = 1./51000. # 1eV-1 = 1/51000*cm
        itg = integrate.quad(lambda theta: self.totXS_intgpart(theta, Ev, A), 0, np.pi)
        tot = 2 * np.pi * GF**2 / (16 * np.pi**2) * (Ev*1e-3)**2 * (self.N - (1-4*0.2325)*self.Z)**2 * itg[0] * 1e-18 * eVminus1tocm**2

        return tot
    

    def diffXS(self, Enr, Ev):
        """
        return differential cross section;
        Ev unit: MeV
        XS unit: cm2/MeV
        """
        GF = 1.1664e-5 # GeV-2
        NA = 6.022e23
        cmtoeVminus1 = 51000 # 1cm = 51000 eV-1
        eVminus1tocm = 1./51000. # 1eV-1 = 1/51000*cm
        xs = GF**2 * (self.mN*1e-9) / (8*np.pi) * ((4*0.2325-1)*self.Z+self.N)**2 * (2 - Enr * self.mN*1e-6/Ev**2) * 1 * 1e-18 * eVminus1tocm**2 * 1e-3
        # firstly omit the form factor
        if xs < 0:
            return 0
        else:
            return xs
