import numpy as np
import matplotlib.pyplot as plt

import scipy.integrate as integrate
from astropy import units as u

from CEvNS_XS import CEvNS_XS
from IBD_XS import IBD_XS
from SNNuloader import SNNuloader
from toyDetector import toyDetector

from tqdm import  tqdm

class visible_spectrum:

    def __init__(self, channel, model, xsec, det) -> None:
        self.channel    = channel
        self.model      = model
        self.xsec       = xsec
        self.det        = det
        self.model.load()

    def getVisibleEventAtEvisAtT(self, t, Evismin, Evismax, mo):

        Evis = (Evismax + Evismin) / 2.
        if self.channel == "IBD":
            deltaE = 1.81 - 1.022
            Evmin, Evmax = 1.81, 80
            Nevt = integrate.quad(lambda Ev: self.det.nH*self.model.getEventAtEarthAtTime(t, Ev, mo, 1)*self.xsec.totXS(Ev), Evmin, Evmax)
            return Nevt[0]


        elif self.channel == "CEvNS":
            Evmin, Evmax, stepEv = 0, 80, 0.1
            NumEv = int((Evmax - Evmin) / stepEv)
            factor0 = 0
            for i in range(NumEv):
                Ev = Evmin + (i+0.5) * stepEv
                for f in range(6):
                    fluence = self.model.getEventAtEarthAtTime(t, Ev, mo, f)
                    factor1 = self.xsec.diffXS(Evis, Ev) * (Evismax - Evismin)
                    #print(t, Evis, Ev, f, fluence, self.xsec.diffXS(Evis, Ev) )
                    factor0 += fluence * factor1 * stepEv
            factor0 = factor0 * self.det.nC 
            return factor0

        elif self.channel == "NuE":
            Evmin, Evmax, stepEv = 0, 80, 0.1
            NumEv = int((Evmax - Evmin) / stepEv)
            factor0 = 0

            for f in range(6):
                val, err = integrate.dblquad(lambda y, x: self.xsec.diffXS(x, y, f)*self.model.getEventAtEarthAtTime(t, y, mo, f), Evismin, Evismax, Evmin, Evmax )
                factor0 += val
            return factor0

            #for i in range(NumEv):
            #    Ev = Evmin + (i+0.5) * stepEv
            #    for f in range(6):
            #        fluence = self.model.getEventAtEarthAtTime(t, Ev, mo, f)
            #        factor1 = self.xsec.diffXS(Evis, Ev, f) * (Evismax - Evismin)
            #        factor0 += fluence * factor1 * stepEv
            #        #print(f, t, Ev, 6*self.det.nC, self.det.nH, fluence, self.xsec.diffXS(Evis, Ev, f), fluence*factor1*stepEv)
            #factor0 = factor0 * (6 * self.det.nC + self.det.nH)
            #return factor0

        else:
            # not implemented yet
            return 0
