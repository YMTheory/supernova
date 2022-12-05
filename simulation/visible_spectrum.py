import numpy as np
import matplotlib.pyplot as plt

import scipy.integrate as integrate
from astropy import units as u
from scipy.misc import derivative
import numba as nb

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
            self.xsec.redefine_ufunc()
        elif channel == "IBD":
            self.xsec   = IBD_XS()
        elif channel == "CEvNS":
            Na = 6.022e23
            gtoeV = 5.61e32
            mC12 = 12 / Na * gtoeV
            self.xsec   = CEvNS_XS(6, 6, mC12)
            self.xsec.redefine_ufunc()

        # Integral of neutrino energy
        self.Evmin = 0
        self.Evmax = 80
        self.stepEv = 0.1
        self.NumEv = int((self.Evmax - self.Evmin) / self.stepEv)


    def getVisibleEventAtTAtEvis(self, t, Evis, mo):

        if self.channel == "NuE":
            npar = self.det.nH + 6 * self.det.nC
            val = 0
            for f in range(6):
                val += integrate.quad(lambda Ev: self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, f) * self.xsec.diffXS_ufunc(Evis, Ev, f), self.Evmin, self.Evmax)[0]
            val = val * npar
            return val

        elif self.channel == "NuP" or channel == "CEvNS":
            if self.channel == "NuP":
                npar = self.det.nH
            else:
                npar = self.det.nC
            val = 0
            for f in range(6):
                val += integrate.quad(lambda Ev: self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, f) * self.xsec.diffXS_ufunc(self.det.proton_unquenchedE(Evis), Ev, f) / derivative(self.det.proton_quenchedE, self.det.proton_unquenchedE(Evis)), self.Evmin, self.Evmax)[0]
            val = val * npar
            return val

        else:
            return 0



    #@nb.jit(nopython=True, parallel=True)
    #@profile
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
                #val += integrate.dblquad(lambda Ev, Evis: self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, f) * self.xsec.diffXS(Evis, Ev, f), Evismin, Evismax, self.Evmin, self.Evmax )[0]
                val += integrate.dblquad(lambda Ev, Evis: npar * self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, f) * self.xsec.diffXS_ufunc(Evis, Ev, f), Evismin, Evismax, self.Evmin, self.Evmax)[0]
            Nevt = val
            return Nevt

        elif self.channel == "NuP" or self.channel == "CEvNS":
            val = 0
            if self.channel == "NuP":
                npar = self.det.nH
            else:
                npar = self.det.nC
            for f in range(6):
                tmp_val, tmp_err= integrate.dblquad(lambda Ev, Evis: self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, f) * self.xsec.diffXS_ufunc(self.det.proton_unquenchedE(Evis), Ev, f) / derivative(self.det.proton_quenchedE, self.det.proton_unquenchedE(Evis), dx=1e-4) , Evismin, Evismax, self.Evmin, self.Evmax)
            Nevt = val * npar
            return Nevt


        else:
            # not implemented yet...
            return 0


    def getVisibleEventAtTEvisIntegral_manually(self, t, Evismin, Evismax, mo):
        """
        do two-fold integration manually
        """
        N_Ev = 80
        stepEv = (self.Evmax - self.Evmin) / N_Ev
        if self.channel == "IBD":
            stepEvis = 1
            npar = self.det.nH
        elif self.channel == "NuE":
            stepEvis = 1
            npar = self.det.nH + 6 * self.det.nC
        elif self.channel == "NuP":
            stepEvis = 0.01
            npar = self.det.nH
        elif self.channel == "CEvNS":
            stepEvis = 0.01
            npar = self.det.nC
        else:
            npar = 0

        N_Evis = int((Evismax - Evismin) / stepEvis)
        
        Nevt = 0

        if self.channel == "IBD":
            deltaE = 1.81 - 1.022
            Evmin_ibd, Evmax_ibd = Evismin + deltaE, Evismax + deltaE
            for iEv in range(N_Ev):
                tmpEv = Evmin_ibd + stepEv * (iEv + 0.5)
                Nevt += npar * self.model.getEventAtEarthAtTime_ufunc(t, tmpEv, mo, 1) * self.xsec.totXS(tmpEv) * stepEv
            return Nevt

        else:
            for iEvis in range(N_Evis):
                tmpEvis = Evismin + stepEvis * (iEvis + 0.5)
                for iEv in range(N_Ev):
                    tmpEv = self.Evmin + stepEv * (iEv + 0.5)

                    if self.channel == "NuE":
                        for f in range(6):
                            val = npar * self.model.getEventAtEarthAtTime_ufunc(t, tmpEv, mo, f) * self.xsec.diffXS_ufunc(tmpEvis, tmpEv, f)
                            #print(t, tmpEvis, tmpEv, f, npar, self.model.getEventAtEarthAtTime_ufunc(t, tmpEv, mo, f), self.xsec.diffXS_ufunc(tmpEvis, tmpEv, f))
                            Nevt += val * stepEv * stepEvis


                    elif self.channel == "NuP" or self.channel == "CEvNS" :
                        for f in range(6):
                            tmpT = self.det.proton_unquenchedE(tmpEvis)
                            dEvis_over_dT = derivative(self.det.proton_quenchedE, tmpT, dx=1e-4)
                            val = npar * self.model.getEventAtEarthAtTime_ufunc(t, tmpEv, mo, f) * self.xsec.diffXS_ufunc(tmpT, tmpEv, f) / dEvis_over_dT
                            Nevt += val * stepEv * stepEvis

                    else:
                        pass
            return Nevt






