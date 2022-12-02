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
        Qw =  self.N - (1-4*0.2325)*self.Z
        print(GF**2/8/np.pi, (Ev*1e-3)**2, Qw**2, itg[0], eVminus1tocm**2)
        tot = 2 * np.pi * GF**2 / (16 * np.pi**2) * (Ev*1e-3)**2 * Qw**2 * itg[0] * 1e-18 * eVminus1tocm**2
        print("totXS = ", tot)

        return tot


    def formFactor1(self, q):
        ### PHYSICS LETTERS B,vol 264, number 1, 2 (1991)
        ### q: unit MeV
        cmtoeVminus1 = 51000 # 1cm = 51000 eV-1
        q = q * 1e6 # eV
        A = self.Z + self.N
        R = 1.2 * np.power(A, 1./3.) * 1e-13 # cm
        s = 1 * 1e-13 # cm
        R0 = np.sqrt(R**2 - 5*s**2)
        
        x = q * R0 * cmtoeVminus1
        j1 = np.sin(x)/x**2 - np.cos(x)/x
        return 3 * j1 / x * np.exp(-0.5 * (q*s*cmtoeVminus1)**2)
        

    
    def diffXS(self, Enr, Ev):
        """
        J. Phys. G, Nucl. Part. Phys. 39, 095204
        return differential cross section;
        Ev unit: MeV
        XS unit: cm2/MeV
        """
        GF = 1.1664e-5 # GeV-2
        NA = 6.022e23
        cmtoeVminus1 = 51000 # 1cm = 51000 eV-1
        eVminus1tocm = 1./51000. # 1eV-1 = 1/51000*cm
        q = np.sqrt(Enr**2 + 2*self.mN*1e-6*Enr)
        xs = GF**2 * (self.mN*1e-9) / (8*np.pi) * ((4*0.2325-1)*self.Z+self.N)**2 * (1 + (1 - Enr / Ev)**2 - self.mN*1e-6*Enr/Ev**2) * self.formFactor1(q)**2 * 1e-18 * eVminus1tocm**2 * 1e-3
        # firstly omit the form factor
        if xs < 0:
            return 0
        else:
            return xs
