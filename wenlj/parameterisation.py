import numpy as np
import matplotlib.pyplot as plt
import ROOT



def getIntegral(hist2d, x1, x2, y1, y2, opt):
    binx1 = hist2d.GetXaxis().FindBin(x1)
    binx2 = hist2d.GetXaxis().FindBin(x2)
    biny1 = hist2d.GetYaxis().FindBin(y1)
    biny2 = hist2d.GetYaxis().FindBin(y2)

    # print('LowEdge(binx1)',  hist2d.GetXaxis().GetBinLowEdge(binx1))
    # print('LowEdge(binx2+1)',  hist2d.GetXaxis().GetBinLowEdge(binx2+1))
    # print('LowEdge(biny1)',  hist2d.GetYaxis().GetBinLowEdge(biny1))
    # print('LowEdge(biny2+1)',  hist2d.GetYaxis().GetBinLowEdge(biny2+1))

    return hist2d.Integral(binx1, binx2, biny1, biny2, opt)



if __name__ == "__main__" :

    modelNum = 82503
    chaname  = 1
    nuType   = -1
    nuMass   = 0.0
    pdfMH    = 1
    dist     = 10.0
    group    = 1
    Ethr     = 0.1


    ## Define Fitting PDF
    fitTmin, fitTmax = -0.03, 0.08
    nuTime = ROOT.RooRealVar("nuTime", "nuTime", fitTmin, fitTmax)
    
    ## 1st exponential component for neutronization burst
    tau1    = ROOT.RooRealVar("tau1", "tau1", 100, 0, 1000)
    shift1  = ROOT.RooRealVar("shift1", "shift1", 1.0, -100, 100)
    nuTime1 = ROOT.RooFormulaVar("nuTime1", "nuTime+shift1", ROOT.RooArgList(nuTime, shift1))

    f1 = ROOT.RooExponential("f1", "1st expo", nuTime1, tau1)

    # 1st step function here 
    nuTimeN = ROOT.RooFormulaVar("nuTimeN", "-nuTime", ROOT.RooArgList(nuTime))
    m0 = ROOT.RooRealVar("m0", "m0", 0.0)
    m1 = ROOT.RooRealVar("m1", "m1", 1.0)
    coefList = ROOT.RooArgList(m0, m1)
    #TArrayD limits(nbins+1)
    limits = ROOT.TArrayD(3)
    limits[0] = -0.08
    limits[1] = -0.01
    limits[2] = 0.04
    sf1 = ROOT.RooParametricStepFunction("step1", "1st step function", nuTimeN, coefList, limits, 2)

    # prod pdf
    prod1 = ROOT.RooProdPdf("prod1", "1st production of expo and step function", ROOT.RooArgSet(f1, sf1))

    ## 2nd exponential component for accretion phase
    tau2    = ROOT.RooRealVar("tau2", "tau2", -100, -1000, 1000)
    shift2  = ROOT.RooRealVar("shift2", "shift2", 0.02, -10, 10)
    nuTime2 = ROOT.RooFormulaVar("nuTime2", "nuTime+shift2", ROOT.RooArgList(nuTime, shift2))

    f2 = ROOT.RooExponential("f2", "2nd expo", nuTime2, tau2)

    # 2nd step function here 
    #TArrayD limits(nbins+1)
    limits[0] = -0.03
    limits[1] = 0.01
    limits[2] = 0.08
    sf2 = ROOT.RooParametricStepFunction("step2", "2nd step function", nuTime, coefList, limits, 2)

    # prod pdf
    prod2 = ROOT.RooProdPdf("prod2", "2nd production of expo and step function", ROOT.RooArgSet(f2, sf2))

    # 3rd unifrom component
    f3 = ROOT.RooUniform("f3", "3rd uniform", nuTime)

    sf3 = ROOT.RooParametricStepFunction("step3", "3rd step function", nuTime, coefList, limits, 2)

    # prod pdf
    prod3 = ROOT.RooProdPdf("prod3", "3rd production of expo and step function", ROOT.RooArgSet(f3, sf3))


    ## Add Three Pdf 
    n1 = ROOT.RooRealVar("n1", "n1", 1.0, 0, 10)
    n2 = ROOT.RooRealVar("n2", "n2", 1.0, 0, 10)
    n3 = ROOT.RooRealVar("n3", "n3", 10.0, 0, 1000)
    model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(prod1, prod2, prod3), ROOT.RooArgList(n1, n2, n3))

    c1 = ROOT.TCanvas()
    xframe = nuTime.frame()
    model.plotOn(xframe)
    xframe.Draw() 
    c1.SaveAs("testPdf.pdf")
