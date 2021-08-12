import matplotlib.pyplot as plt
import numpy as np
import ROOT
from array import array
import namespace as ns

libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)

imode = 82503
dist = 10
nuType = 0
MH = 2
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

EvisTmp = 12
tTmp = 0.03

def EvisSpectrumAtTime():
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
            sum_weight += Npar*fluence*dxs*evisFactor    # simple sum for integral

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


   
def EvisSpectrum(tTmp,  iMH):
    Evmin, Evmax, step_Ev = 0.0, 80.0, 0.5
    nbins_Ev = round((Evmax - Evmin)/step_Ev)

    Evismin, Evismax, step_Evis = 0.1, 25, 0.02
    nbins_Evis = round((Evismax - Evismin)/step_Evis)
    
    timeTmp = array('d', [-999.]) # unit: s

    # a specified instance
    timeTmp[0] = tTmp
    
    flux_Evis = []
    for i in range(nbins_Evis):
        EvisTmp = Evismin + step_Evis * (i+0.5)
        me = 0.511
        EvminTmp = 0.5*(EvisTmp + np.sqrt(EvisTmp**2 + 2*EvisTmp*me))
        EvmaxTmp = modelSrc.getSNmaxEv()

        sum_weight = 0.
        for j in range(nbins_Ev):
            EvTmp = Evmin + step_Ev*(j+0.5)
            if EvTmp < EvminTmp:
                continue
            fluence = modelSrc.snFluenceDetAtTime(timeTmp, nuMass, EvTmp, nuType, iMH)
            dxs = peffLS.differentialXS(EvTmp, EvisTmp, 0)

            if fluence>0. and dxs>0. :
                #print("weight from generatePDFs: %.3f" %(Npar*fluence*dxs))
                evisFactor = step_Ev 
                #evisFactor = 1
                sum_weight += Npar*fluence*dxs*evisFactor    # simple sum for integral

        flux_Evis.append(sum_weight)

    sum_flux = 0
    for i in flux_Evis:
        sum_flux += i * step_Evis


    print("Time integral flux from generatePDFs at this time %.3f s : %.3f" %(timeTmp[0], sum_flux) )
    #print("Time integral flux from SNdetect at this time : %.3f" %(snDet.getEventAboveEthrVisAtTime(tTmp, Ethr, nuType, MH)))

    return sum_flux


def loadPDFfile(filename, histname):
    fin = ROOT.TFile(filename, "read")
    hh  = fin.Get(histname)

    # projection manually
    cent_T, cont_T = [], []
    for i in range(hh.GetNbinsX()):
        cent_T.append(hh.GetXaxis().GetBinCenter(i+1))
        yy = 0
        for j in range(hh.GetNbinsY()):
            yy += hh.GetBinContent(i+1, j+1) * hh.GetYaxis().GetBinWidth(j+1)

        cont_T.append(yy)

    return cent_T, cont_T




def TimeSpectra():

    mod = 82503
    dist = 10
    nuMass = 0.0
    Evismax = 25.00
    chaname = 1
    nuType = 0

    #tt = np.arange(0.05, 0.1, 0.001)
    #arr0, arr1, arr2 = [], [], []
    #for i in tt:
    #    tTmp = i
    #    arr0.append(EvisSpectrum(tTmp, 0))
    #    arr1.append(EvisSpectrum(tTmp, 1))
    #    arr2.append(EvisSpectrum(tTmp, 2))
    #    #arr2.append(snDet.getEventAboveEthrVisAtTime(tTmp, Ethr, nuType, MH))
    #    print(i, arr0[-1], arr1[-1], arr2[-1])

    #cent_T, cont_T = loadPDFfile("/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/TEvisPDF_mod82503_cha1nue_nuType0_mh0_mNu0.0eV_10.0kpc_0.1s_Evismax25.root", "hET_mod82503_cha1_mh0")
    tt, arr0 = loadPDFfile(ns.pdfSumFileName(mod, dist, chaname ,nuMass, 0, nuType, Evismax), "TEvisPdf")
    tt, arr1 = loadPDFfile(ns.pdfSumFileName(mod, dist, chaname ,nuMass, 1, nuType, Evismax), "TEvisPdf")
    tt, arr2 = loadPDFfile(ns.pdfSumFileName(mod, dist, chaname ,nuMass, 2, nuType, Evismax), "TEvisPdf")
    
    plt.plot(tt, arr0, "-" , label="no osc", color="royalblue")
    plt.plot(tt, arr1, "-" , label="NO", color="orange")
    plt.plot(tt, arr2, "-" , label="IO", color="seagreen")


    nuType = -1
    tt, arr3 = loadPDFfile(ns.pdfSumFileName(mod, dist, chaname ,nuMass, 0, nuType, Evismax), "TEvisPdf")
    tt, arr4 = loadPDFfile(ns.pdfSumFileName(mod, dist, chaname ,nuMass, 1, nuType, Evismax), "TEvisPdf")
    tt, arr5 = loadPDFfile(ns.pdfSumFileName(mod, dist, chaname ,nuMass, 2, nuType, Evismax), "TEvisPdf")
    
    plt.plot(tt, arr3, "--" , label="no osc", color="royalblue")
    plt.plot(tt, arr4, "--" , label="NO", color="orange")
    plt.plot(tt, arr5, "--" , label="IO", color="seagreen")


    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"dN/dt($s^{-1}$)")
    plt.savefig("./spectra/valid/eESRate_osc_10kpc_0.0eV_onlynue+allflavour.pdf")
    plt.show()


def integral(hh, tlow, thigh, Elow, Ehigh):
    xlow = hh.GetXaxis().FindBin(tlow)
    xhig = hh.GetXaxis().FindBin(thigh)
    ylow = hh.GetYaxis().FindBin(Elow)
    yhig = hh.GetYaxis().FindBin(Ehigh)

    return hh.Integral(xlow, xhig, ylow, yhig, "width")





if __name__ == "__main__":
    #EvisSpectrumAtTime()
    #numIntg()
    #EvisSpectrum()

    TimeSpectra()
    #print(snDet.getEventAboveEthrVis(0.1, 0, 0))
    #print(snDet.getEventAboveEthrVis(0.1, 0, 1))
    #print(snDet.getEventAboveEthrVis(0.1, 0, 2))
    #print(snDet.getEventAboveEthrVisTimeInterval(-0.1, 0.1, 0.1, 0, 0))
    #print(snDet.getEventAboveEthrVisTimeInterval(-0.1, 0.1, 0.1, 0, 1))
    #print(snDet.getEventAboveEthrVisTimeInterval(-0.1, 0.1, 0.1, 0, 2))


    #f1 = ROOT.TFile("/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/TEvisPDF_mod82503_cha1nue_mh0_mNu0.0eV_10.0kpc_0.1s_Evmax25.root", "read")
    #h1 = f1.Get("hET_mod82503_cha1_mh0")
    #f2 = ROOT.TFile("/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/TEvisPDF_mod82503_cha1nue_nuType0_mh0_mNu0.0eV_10.0kpc_0.1s_Evismax25.root", "read")
    #h2 = f2.Get("hET_mod82503_cha1_mh0")
    #print("NEvents from f1:" , integral(h1, -0.04, 0.1, 0.1, 10))
    #print("NEvents from f2" , integral(h2, -0.04, 0.1, 0.1, 10))





