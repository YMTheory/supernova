import numpy as np
import matplotlib.pyplot as plt

import scipy.integrate as integrate
from astropy import units as u

from CEvNS_XS import CEvNS_XS
from IBD_XS import IBD_XS
from NuE_XS import NuE_XS
from NuP_XS import NuP_XS
from SNNuloader import SNNuloader
from toyDetector import toyDetector

from tqdm import  tqdm

class visible_spectrum:

    def __init__(self, channel, SNmass=25, EoS="LS220", dist=10, detMass=20000, fracH=0.12, fracC=0.88) -> None:
        self.channel    = channel

        self.model      = SNNuloader(SNmass, "Accretion", EoS, dist)
        self.model.load()
        self.model.redefine_ufunc()

        self.det        = toyDetector(detMass, fracH, fracC)
        
        if channel == "NuE":
            self.xsec   = NuE_XS()
            self.xsec.redefine_ufunc()
        elif channel == "NuP":
            self.xsec   = NuP_XS()
        elif channel == "IBD":
            self.xsec   = IBD_XS()
        elif channel == "CEvNS":
            Na = 6.022e23
            gtoeV = 5.61e32
            mC12 = 12 / Na * gtoeV
            self.xsec   = CEvNS_XS(6, 6, mC12)

        # Integral of neutrino energy
        self.Evmin = 0
        self.Evmax = 80
        self.stepEv = 0.1
        self.NumEv = int((self.Evmax - self.Evmin) / self.stepEv)


    def getVisibleEventAtTEvisIntegral(self, t, Evismin, Evismax, mo):

        if self.channel == "IBD":
            deltaE = 1.81 - 1.022
            Enumin = Evismin + deltaE
            Enumax = Evismax + deltaE
            Nevt = integrate.quad(lambda Ev: self.det.nH * self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, 1) * self.xsec.totXS(Ev), Enumin, Enumax )[0]
            return Nevt

        elif self.channel == "NuE":
            val = 0
            npar = self.det.nH + 6 * self.det.nC
            for f in range(6):
                val += integrate.dblquad(lambda Ev, Evis: self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, f) * self.xsec.diffXS(Evis, Ev, f), Evismin, Evismax, self.Evmin, self.Evmax )[0]
            Nevt = npar * val
            return Nevt

        elif self.channel == "NuP" or self.channel == "CEvNS":
            val = 0
            if self.channel == "NuP":
                npar = self.det.nH
            else:
                npar = self.det.nC
            Tmin = self.det.proton_unquenchedE(Evismin)
            Tmax = self.det.proton_unquenchedE(Evismax)
            #print(f"[{Evismin}, {Evismax}] -> [{Tmin}, {Tmax}]")
            for f in range(6):
                val += integrate.dblquad(lambda Evis, Ev: self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, f) * self.xsec.diffXS(self.det.proton_unquenchedE(Evis), Ev, f), self.Evmin, self.Evmax, Evismin, Evismax)[0]
            Nevt = npar * val
            return Nevt


        else:
            # not implemented yet...
            return 0



        





