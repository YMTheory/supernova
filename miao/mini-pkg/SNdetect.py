import numpy as np
import SNGarchingIntegFcn as gar
from SNnumGarchingSrc import SNnumGarchingSrc
import SNnueXS as nuexs
from detector import detector
import ROOT

class SNdetect(object):
    
    mode   = 82503
    dist   = 10
    fEvmin = 0
    fEvmax = 60
    SNGar  = SNnumGarchingSrc(mode, dist) 
    det    = detector()



    def __init__(self, imode, d):
        self.mode = imode
        self.dist = d

        self.SNGar.setDistance(d)
        self.SNGar.setModeNum(d)


    """
    def getEventAboveEthrVisAtTime(time, Ethr, tp, MH):
        fEvismax = 0
        ## below is for eES channel
        me = 0.511
        fEvismax = fEvmax/(1+me/(2.*fEvmax))

        f = ROOT.TF1("eveVisatTime", fcnEventsVisatTime, 0, fEvismax, 3)
        f.SetParameter(0, time)
        f.SetParameter(1, tp)
        f.SetParameter(2, MH)
    
    """    





    # -------- FCN -------- #
    def fcnAnaEv2TatTime(self, x, p):
        Ev = x[0]
        time = p[0]
        T = p[1]
        tp = p[2]
        MH = p[3]

        ntype = 6
        totflu = 0
        if tp == -1:
            for it in range(ntype):
                flu = SNGar.oneSNFluenceDetTimeShift(time, Ev, 0, it, MH)
                dxs = nuexs.differentialXS(Ev, T, it)
                totflu += flu * dxs 
        else:
            flu = SNGar.oneSNFluenceDetMH(time, Ev, 0, ty, MH)
            dxs = nuexs.differentialXS(Ev, T, ty)
            totflu += flu * dxs


        return totflu


    def getTSpectrumAtTime(self, time, T, tp, MH):
        Evmin = 0
        Npar = 0
        me = 0.511
        Evmin = 0.5 * np.sqrt(T+np.sqrt(T**2+2*T*me))
        Npar = self.det.getNumberOfProton() + 6*self.det.getNumberOfCarbon()

        f = ROOT.TF1("anaEv2TatTime", self.fcnAnaEv2TatTime, self.fEvmin, self.fEvmax, 4)
        f.SetParameter(0, time)
        f.SetParameter(1, T)
        f.SetParameter(2, tp)
        f.SetParameter(3, MH)
        wf1 = ROOT.Math.WrappedTF1(f)
        ig = ROOT.Math.GSLIntegrator(ROOT.Math.IntegrationOneDim.kADAPTIVE,ROOT.Math.Integration.kGAUSS61);
        ig.SetFunction(wf1);
        ig.SetRelTolerance(1e-3);
        spect = 0
        ig.Integral()
        #spect = Npar * ig.Integral(Evmin, self.fEvmax)

        return spect


    def getTFromEvis(self, Evis):
        return Evis


    def getEvisSpectrumAtTime(self, time, Evis, tp, MH):
        return self.getTSpectrumAtTime(time, self.getTFromEvis(Evis), tp, MH)



    def fcnEventsVisatTime(self, x, par):
        Evis = x[0]
        time = par[0]
        tp = int(par[1])
        MH = int(par[2])

        return self.getEvisSpectrumAtTime(time, Evis, tp, MH)



    ## analytical energy distribution
    def getEventAboveEthrVisAtTime(self, time, Ethr, tp, MH):
        Evismax = 0
        # only for eES channel here ...
        me = 0.511
        Evismax = self.fEvmax / (1+me/(2.*self.fEvmax))
        
        f = ROOT.TF1("eveVisatTime", fcnEventsVisatTime, 0, Evismax, 3)




















