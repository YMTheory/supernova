import matplotlib.pyplot as plt
import numpy as np
import ROOT
from array import array

libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)

imode = 82503
dist = 10
nuType = 0
MH = 1
nuMass = 0.0

# === Configuration === #
snDet = ROOT.SNdetect.instance()
snDet.setSrcModel(imode)
modelSrc = snDet.getPointerSrc()

snDet.initChannel(1)
modelSrc.setSNDistance(dist)
Ethr = 0.1
snDet.getPointerEffectLS().setThresholdE(Ethr)
snDet.initFCN()
peffLS = snDet.getPointerEffectLS()
Npar = peffLS.getNumberOfProton()+6*peffLS.getNumberOfCarbon()

# === === === === === === #

EvisTmp = 10
tTmp = 0.02

def EvisSpectrum():
    Evmin, Evmax, step_Ev = 0.0, 80.0, 0.5
    nbins_Ev = round((Evmax - Evmin)/step_Ev)

    Evismin, Evismax, step_Evis = 0.1, 25, 0.02
    nbins_Evis = round((Evismax - Evismin)/step_Evis)
    
    timeTmp = array('d', [-999.]) # unit: s

    # a specified instance
    timeTmp[0] = tTmp
    
    me = 0.511
    EvminTmp = 0.5*(EvisTmp + np.sqrt(EvisTmp**2 + 2*EvisTmp*me))
    EvmaxTmp = modelSrc.getSNmaxEv()

    sum_weight = 0.
    for j in range(nbins_Ev):
        EvTmp = Evmin + step_Ev*(j+0.5)
        if EvTmp < EvminTmp:
            continue
        fluence = modelSrc.snFluenceDetAtTime(timeTmp, nuMass, EvTmp, nuType, MH)
        dxs = peffLS.differentialXS(EvTmp, EvisTmp, 0)

        if fluence>0. and dxs>0. :
            #print("weight from generatePDFs: %.3f" %(Npar*fluence*dxs))
            evisFactor = step_Ev 
            #evisFactor = 1
            sum_weight += Npar*fluence*dxs*evisFactor

    print("Total weight from generatePDFs : %.3f" %(sum_weight))

    print("Total weight from SNdetect  %.3f" %(snDet.getEvisSpectrumAtTime(timeTmp[0], EvisTmp, nuType, MH)))


def numIntgFCN(Ev):

    time = tTmp
    T = EvisTmp

    flu = snDet.getPointerSrc().oneSNFluenceDetAtTime(time, Ev, nuType, MH)
    dxs = snDet.getPointerEffectLS().differentialXS(Ev, T, nuType)
    totflu = flu * dxs * Npar

    return totflu


from scipy import integrate
def numIntg():

    me = 0.51
    EvminTmp = 0.5*(EvisTmp + np.sqrt(EvisTmp**2 + 2*EvisTmp*me))
    EvmaxTmp = modelSrc.getSNmaxEv()
    
    A, e = integrate.quad(numIntgFCN, EvminTmp, EvmaxTmp)
    print(A, e)

    



if __name__ == "__main__":
    EvisSpectrum()
    numIntg()
