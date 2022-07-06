import numpy as np
import uproot as up
from iminuit import Minuit
import sys
import ROOT
import math
import matplotlib.pyplot as plt
plt.style.use("science")

def getIntegral(hist1d, x1, x2, opt):
    binx1 = hist1d.GetXaxis().FindBin(x1)
    binx2 = hist1d.GetXaxis().FindBin(x2)
    return hist1d.Integral(binx1, binx2, opt)

def getIntegral(hist1d, x1, x2, opt):
    binx1 = hist1d.GetXaxis().FindBin(x1)
    binx2 = hist1d.GetXaxis().FindBin(x2)
    return hist1d.Integral(binx1, binx2, opt)


def calc_NLL(data, pdf, dT)->float:
    nll = 0
    for i in data:
        nll += np.log(pdf.Interpolate(i+dT))
    nu = getIntegral(pdf, -20+dT, 20+dT, "width")  # expected event numebr
    N = len(data)                                  # observed event numebr

    #nll += N * np.log(nu)
    #nll -= nu

    return -nll



if __name__ == "__main__" :

    mod = 82503
    MH = "IO"
    cha = "pES"
    Ethr = 0.2

    if len(sys.argv) > 1:
        ### parametere configuration
        print("argument number: %d"%int((len(sys.argv)-1)/2))
        for i in range(1, len(sys.argv)-1, 2):
            print("==> Fitting ")
            if sys.argv[i] == "-model":
                mod = int(sys.argv[i+1])
                print("=====> Model: %d"%mod)
            elif sys.argv[i] == "-cha" :
                cha = sys.argv[i+1]
                print("=====> Channel: %s"%cha)
            elif sys.argv[i] == "-mo" :
                MH = sys.argv[i+1]
                print("=====> Mass ordering: %s"%MH)
            elif sys.argv[i] == "-Ethr" :
                Ethr = float(sys.argv[i+1])
                print("=====> Ethr: %.2f MeV"%Ethr)
            else:
                print("Error: No such argument !")
                exit(-1)


    ## read pdf
    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Garching%d_PDF_NO_10kpc_%s_%.2fMeV.root" %(mod, cha, Ethr)
    hn = cha
    pdff1 = ROOT.TFile(pdffile, "read")
    pdf1 = pdff1.Get(hn)
    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Garching%d_PDF_IO_10kpc_%s_%.2fMeV.root" %(mod, cha, Ethr)
    hn = cha
    pdff2 = ROOT.TFile(pdffile, "read")
    pdf2 = pdff2.Get(hn)

    Tmin, Tmax = -20, 20
    nuTime1D = ROOT.RooRealVar("nuTime1D", "nuTime1D", Tmin, Tmax)
    nuTime1D.setRange('fullRange', Tmin, Tmax)
    sigHist1 = ROOT.RooDataHist('sigHist1', 'sigHist1', ROOT.RooArgList(nuTime1D), ROOT.RooFit.Import(pdf1, True))
    sigHist2 = ROOT.RooDataHist('sigHist2', 'sigHist2', ROOT.RooArgList(nuTime1D), ROOT.RooFit.Import(pdf2, True))
    sigPdf1 = ROOT.RooHistPdf('sigPdf1', 'sigPdf1', ROOT.RooArgList(nuTime1D), sigHist1, 2)
    sigPdf2 = ROOT.RooHistPdf('sigPdf2', 'sigPdf2', ROOT.RooArgList(nuTime1D), sigHist2, 2)
    NEVENTS1 = getIntegral(pdf1, -20, 20, "width")
    NEVENTS2 = getIntegral(pdf2, -20, 20, "width")
    print("NEVENTS (NO): %d, NEVENTS (IO): %d" %(NEVENTS1, NEVENTS2))
    nsig1 = ROOT.RooRealVar("nsig1","#signal events", NEVENTS1, 0., NEVENTS1*10000)
    nsig2 = ROOT.RooRealVar("nsig2","#signal events", NEVENTS2, 0., NEVENTS2*10000)
    fitPdf1  = ROOT.RooExtendPdf('fitPdf1', 'fitPdf1', sigPdf1, nsig1)
    fitPdf2  = ROOT.RooExtendPdf('fitPdf2', 'fitPdf2', sigPdf2, nsig2)
    
    # declare variables
    Tmin, Tmax = -20, 20     # Fitting time range , unit : ms
    datafile = "/junofs/users/miaoyu/supernova/production/Data/10kpc/Garching%d_%s_data_%s_10kpc_thr%.2fMeV.root" %(mod, cha, MH, Ethr)
    print(datafile)
    dataFile = ROOT.TFile(datafile, "red")
    inputTree = dataFile.Get('tFixedStat')

    nuTime = ROOT.RooRealVar("nuTime", "nuTime", Tmin, Tmax)
    evtID    = ROOT.RooRealVar("evtID", "evtID", 0, 500)

    newArgSet = ROOT.RooArgSet(nuTime)

    nStat = 5
    for iSub in range(nStat):
        # Get RooDataSet from the tree
        preSelData1D = ROOT.RooDataSet("dataGroup%d"%(iSub), "dataset with (e,t)", 
                        inputTree, ROOT.RooArgSet(nuTime, evtID), 
                        "evtID==%d"%(iSub) )

        nFitEvts1D = round(preSelData1D.sumEntries())
        print('preSelData1D.sumEntries(): ', nFitEvts1D )

        ##### Likelihood Coarse scan ##### 
        nllList1, tScan1, nsigList1 = [], [], []
        nllList2, tScan2, nsigList2 = [], [], []
        extDTmax, extDTmin = 10, -10
        locMinNLL1, bestfitDT1 = 100.0, 1.0
        locMinNLL2, bestfitDT2 = 100.0, 1.0
        for dT in range(-3, 9, 1):
            if dT > extDTmax or dT < extDTmin:
                print("Warning! : data might be shifted outside the fitting range !!!")
            fitdata = ROOT.RooDataSet('fitdata', 'fitdata', newArgSet)
            nuTime1D.setRange('nllRange', Tmin+dT, Tmax+dT)
            nuTime.setRange('nllRange', Tmin+dT, Tmax+dT)

            for iEvt in range(nFitEvts1D):
                argSet      = preSelData1D.get(iEvt)
                timeTmp     = argSet.getRealValue('nuTime')
                weight      = preSelData1D.weight()
                weightError = preSelData1D.weightError()

                newArgSet.setRealValue('nuTime1D', timeTmp + dT)

                fitdata.add( newArgSet, weight, weightError)
    
            nuTime1D.setRange('nllRange', Tmin+dT, Tmax+dT)
            res = fitPdf1.fitTo(fitdata, ROOT.RooFit.Range('nllRange'), 
                                ROOT.RooFit.Save(),
                                ROOT.RooFit.Extended(),
                                ROOT.RooFit.PrintLevel(-1))
            nllVal = res.minNll()
            
            if math.isnan(nllVal) != True:
                tScan1.append(dT)
                nllList1.append( nllVal )
                nsigList1.append(nsig1.getVal())
                if locMinNLL1 > nllVal:
                    locMinNLL1 = nllVal
                    bestfitDT1 = dT
            res = fitPdf2.fitTo(fitdata, ROOT.RooFit.Range('nllRange'), 
                                ROOT.RooFit.Save(),
                                ROOT.RooFit.Extended(),
                                ROOT.RooFit.PrintLevel(-1))
            nllVal = res.minNll()
            del fitdata
            
            print("Fitting status ", res.status())
            if math.isnan(nllVal) != True:
                tScan2.append(dT)
                nllList2.append( nllVal )
                nsigList2.append(nsig1.getVal())
                if locMinNLL2 > nllVal:
                    locMinNLL2 = nllVal
                    bestfitDT2 = dT
                
        
        print(tScan1)
        print(nllList1)
        print(nsigList1)
        print(tScan2)
        print(nllList2)
        print(nsigList2)

