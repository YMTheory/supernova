import numpy as np
import ROOT

from CEvNS_XS import CEvNS_XS
from IBD_XS import IBD_XS
from NuE_XS import NuE_XS
from NuP_XS import NuP_XS
from SNNuloader import SNNuloader
from toyDetector import toyDetector

from tqdm import  tqdm


class integration:

    def __init__(self, channel, model, det) -> None:
        print("****** A object of integration has been constructed here!")
        self.channel    = channel

        self.model      = model
        self.model.load()
        self.model.redefine_ufunc()

        self.det        = det
        
        if channel == "NuE":
            self.xsec   = NuE_XS()
            self.xsec.redefine_ufunc()
            self.npar = self.det.nH + 6 * self.det.nC
        elif channel == "NuP":
            self.xsec   = NuP_XS()
            self.xsec.redefine_ufunc()
            self.npar = self.det.nH 
        elif channel == "IBD":
            self.xsec   = IBD_XS()
            self.npar = self.det.nH 
        elif channel == "CEvNS":
            Na = 6.022e23
            gtoeV = 5.61e32
            mC12 = 12 / Na * gtoeV
            self.xsec   = CEvNS_XS(6, 6, mC12)
            self.xsec.redefine_ufunc()
            self.npar = self.det.nC

        # Integral of neutrino energy
        self.Evmin = 0
        self.Evmax = 80



    def __call__(self, x, par):
        Evis     = x[0]
        Ev       = x[1]
        t        = par[0]
        moID     = par[1]
        flavorID = par[2]
        if moID == 0: # NO
            mo = "no"
        elif moID == 1: # ID
            mo = "io"
        return self.npar * self.model.getEventAtEarthAtTime_ufunc(t, Ev, mo, flavorID) * self.xsec.diffXS_ufunc(Evis, Ev, flavorID)



    def __del__(self):
        print("****** Object has been deleted!")



