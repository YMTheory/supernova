import numpy as np
import matplotlib.pyplot as plt
import SNGarchingIntegFcn as gar

import ROOT

class SNnumGarchingSrc(object):

    dist = 10   #kpc
    mode = 85200
    loadFlag = False
    datapath = "/junofs/users/miaoyu/supernova/wenlj/simulation/data/Garching/"
    grspectrum = []

    def __init__(self, imode, d):
        self.dist = d
        self.mode = imode
        self.loadFlag = False
        gar.readFluxGraph(self.mode)

    def readSpectrum(self):
        prefix = "energySpeFiles/energySpecGarching_"
        fullname = self.datapath+prefix+str(self.mode)+".root"
        fin = ROOT.TFile(fullname, "read")
        if not fin:
            print("No energy spectrum file %s !" %fullname)

        for i in range(3):
            g = fin.Get("grmod"+str(self.mode)+"type"+str(i))
            self.grspectrum.append(g)
        
        self.loadFlag = True


    def setDistance(self, d):
        self.dist = d


    def oneSNFluenceDet(self, E, tp):
        if self.loadFlag == False:
            self.readSpectrum()
        index = 1/(4*np.pi)/3.086e21**2/self.dist/self.dist
        if tp < 2:
            return index * self.grspectrum[tp].Eval(E)
        if tp >= 2:
            return index * self.grspectrum[2].Eval(E)
        return 0



    def oneSNFluenceDetMH(self, E, tp, MH):
        # consider oscillation
        fluence = 0
        if MH == 0:
            return self.oneSNFluenceDet(E, tp)
        p, pbar = 0, 0
        if MH == 1:
            p = 0.022
            pbar = 0.687
        if MH == 2:
            p = 0.291
            pbar = 0.022
        if tp == 0:
            fluence = p*self.oneSNFluenceDet(E, 0) + (1-p)*self.oneSNFluenceDet(E, 2)   # nu_e
        if tp == 1:
            fluence = pbar*self.oneSNFluenceDet(E, 1) + (1 - pbar) * self.oneSNFluenceDet(E, 3)
        if tp == 2 or tp == 4:
            fluence = 0.5*(1-p)*self.oneSNFluenceDet(E, 0) + 0.5*(1+p)*self.oneSNFluenceDet(E, 2)
        if tp ==3 or tp == 5:
            fluence = 0.5*(1-pbar)*self.oneSNFluenceDet(E, 0) + 0.5*(1+pbar)*self.oneSNFluenceDet(E, 2)

        return fluence



    def totalSNFluenceDet(self, E, MH):
        tot_tp = 6
        tfluence = 0
        for i in range(tot_tp):
            tfluence += self.oneSNFluenceDetMH(E, i, MH)
        return tfluence
            

    def oneSNFluenceDetTimeIntegE(self, time, tp, MH):
        index = 1/(4*np.pi)/3.086e21**2/self.dist/self.dist
        fluence = 0
        if MH == 0:
            return gar.getNumT(time, tp)
        p, pbar = 0, 0
        if MH == 1:
            p = 0.022
            pbar = 0.687
        if MH == 2:
            p = 0.291
            pbar = 0.022
        if tp == 0:
            fluence = p*gar.getNumT(time, 0) + (1-p)*gar.getNumT(time, 2)   # nu_e
        if tp == 1:
            fluence = pbar*gar.getNumT(time, 1) + (1 - pbar) * gar.getNumT(time, 3)
        if tp == 2 or tp == 4:
            fluence = 0.5*(1-p)*gar.getNumT(time, 0) + 0.5*(1+p)*gar.getNumT(time, 2)
        if tp ==3 or tp == 5:
            fluence = 0.5*(1-pbar)*gar.getNumT(time, 0) + 0.5*(1+pbar)*gar.getNumT(time, 2)

        return fluence*index
        


    def totalSNFluenceDetTimeIntegE(self, time, MH):
        tot_tp = 6
        tfluence = 0
        for i in range(tot_tp):
            tfluence += self.oneSNFluenceDetTimeIntegE(time, i, MH)
        return tfluence


    # MODIFY TIME PROFILE DUE TO NEUTRINO MASSES
    def oneSNFluenceDetTimeShift(self, time, nuMass, E, tp, MH):
        index = 1/(4*np.pi)/3.086e21**2/self.dist/self.dist
        fluence = 0
        if MH == 0 :
            fluence = gar.getEventAtTime(time, E, tp)
        p, pbar = 0, 0
        if MH == 1:
            p = 0.022
            pbar = 0.687
        if MH == 2:
            p = 0.291
            pbar = 0.022
        if tp == 0:
            fluence = p*gar.getEventAtTime(time, E, 0) + (1-p)*gar.getEventAtTime(time, E, 2)   # nu_e
        if tp == 1:
            fluence = pbar*gar.getEventAtTime(time, E, 1) + (1 - pbar) * gar.getEventAtTime(time,  E, 3)
        if tp == 2 or tp == 4:
            fluence = 0.5*(1-p)*gar.getEventAtTime(time, E, 0) + 0.5*(1+p)*gar.getEventAtTime(time, E, 2)
        if tp ==3 or tp == 5:
            fluence = 0.5*(1-pbar)*gar.getEventAtTime(time, E, 0) + 0.5*(1+pbar)*gar.getEventAtTime(time, E, 2)

        # modification :
        DeltaT = 5.14e-3 * nuMass**2 * 100./E**2 * (self.dist/10)     # unit: s

        return fluence*index, time+DeltaT



    def totalSNFluenceDetTimeShift(self, time, nuMass, E, MH):
        tot_tp = 6
        tflux = 0
        for i in range(tot_tp):
            flux, timeCorr = self.oneSNFluenceDetTimeShift(time, nuMass, E, i, MH)
            tflux += flux
        return tflux, timeCorr






