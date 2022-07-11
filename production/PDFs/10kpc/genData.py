import numpy as np
import matplotlib.pyplot as plt
import ROOT
from array import array


def loadPDF(mod, cha, MO, thr):
    filename = "Garching%d_timespec_%s_10kpc_thr%.1fMeV.root" %(mod, MO, thr)
    f = ROOT.TFile(filename, "read")
    h = f.Get(cha)

    x, y = [], []
    N = h.GetNbinsX()
    for i in range(N):
        x.append(h.GetBinCenter(i+1))
        y.append(h.GetBinContent(i+1))
    x = np.array(x)
    y = np.array(y)

    return x, y

def getIntegral(hist1d, x1, x2, opt):
    binx1 = hist1d.GetXaxis().FindBin(x1)
    binx2 = hist1d.GetXaxis().FindBin(x2)
    return hist1d.Integral(binx1, binx2, opt)


if __name__ == "__main__" :


    #for imod in [81120,81121,81122,81123,82500,82501,82502,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502,92503,92700,92701,92702,92703,94000,94001,94002,94003]:
    #for imod in [6500]:
    for imod in [82503]:
        for imh in [1, 2]:
            for cha in [0, 1, 2]:
            #for cha in [0]:
                if cha == 0 :
                    Ethr = [0.10, 0.15, 0.20]
                if cha == 1 or cha ==2 :
                    Ethr = [0.20]
                if imh == 1:
                    MO = "NO"
                elif imh == 2:
                    MO = "IO"
                chaname = ["pES", "eES", "IBD"]

                for E in Ethr:
                    filename = "Garching82503_PDF_%s_10kpc_%s_%.2fMeV.root"%( MO, chaname[cha], E)
                    #filename = "Burrows2D12_PDF_%s_10kpc_%s_%.2fMeV.root"%( MO, chaname[cha], E)
                    dataFile = ROOT.TFile(filename, "read")
                    inputHist1D = dataFile.Get("h1")


                    Tmin = inputHist1D.GetXaxis().GetXmin()
                    Tmax = inputHist1D.GetXaxis().GetXmax()
                    print("Tmin = %.1f ms, Tmax = %.1f ms"%(Tmin, Tmax))
                    nuTime   = ROOT.RooRealVar("nuTime", "nuTime", Tmin, Tmax)
                    nuTime.setRange('fullRange', Tmin, Tmax)
                    NEVENTS = getIntegral(inputHist1D, Tmin, Tmax, "width")
                    print("NEVENTS in fullRange : %.2f" %NEVENTS)

                    Tmin, Tmax = -20, 20
                    subNEVENTS = getIntegral(inputHist1D, Tmin, Tmax, 'width')
                    print('NEVENTS in ROI : ', subNEVENTS)
                    #subNEVENTS = round(subNEVENTS)
                    subNEVENTS = int(round(subNEVENTS))
                    print('NEVENTS in ROI (round): ', subNEVENTS)
                    nuTime.setRange('preSelRange', Tmin, Tmax)
                    
                    dataHist1D = ROOT.RooDataHist('dataHist1D', 'dataHist1D',
                            ROOT.RooArgList(nuTime),
                            ROOT.RooFit.Import(inputHist1D, True))
                    dataPdf1D = ROOT.RooHistPdf("dataPdf", "dataPdf",
                            ROOT.RooArgSet(nuTime), dataHist1D, 2)
                    

                    filename = "../../Data/10kpc/Garching82503_%s_data_%s_10kpc_thr%.2fMeV.root" %( chaname[cha], MO, E)
                    print(filename)
                    outFile = ROOT.TFile(filename, "recreate")
                    evtID = array('i', [0])
                    time = array('d', [0.])
                    tFixedStat = ROOT.TTree("tFixedStat", "data set")
                    tFixedStat.Branch("evtID", evtID, "evtID/I")
                    tFixedStat.Branch("nuTime1D", time, 'nuTime/D')

                    sampleSize = 500
                    for i in range(sampleSize):
                        #if i%10 == 0:
                            #print("Running event %d ..."%i)
                        evtID[0] = i
                        data1D = dataPdf1D.generate(ROOT.RooArgSet(nuTime), 3*NEVENTS)
                        preSelData1D = data1D.reduce(ROOT.RooFit.CutRange('preSelRange'))

                        for k in range(subNEVENTS):
                            argSet  = preSelData1D.get(k)
                            time[0] =  argSet.getRealValue('nuTime')

                            tFixedStat.Fill()

                        del data1D
                        del preSelData1D

                    tFixedStat.Write()

                    outFile.Close()


