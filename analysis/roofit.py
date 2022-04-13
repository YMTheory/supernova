import sys
import ROOT
import math
import uproot3
import uproot as up

def loadInit(mod, mh):
    filename = "/junofs/users/miaoyu/supernova/analysis/initValue_2kpc100_mh%d.root"%(mh)
    print(filename)
    ff = up.open(filename)
    model = ff["fit"]["model"].array()
    p0    = ff["fit"]["p0"].array() 
    p1    = ff["fit"]["p1"].array() 
    p2    = ff["fit"]["p2"].array() 
    p3    = ff["fit"]["p3"].array() 
    p4    = ff["fit"]["p4"].array() 

    for i, j, k, x, y, z in zip(model, p0, p1, p2, p3, p4):
        if i == mod:
            return j, k, x, y, z
    print("Wrong model ID!")
    return -1 , -1, -1, -1, -1


if __name__ == "__main__" :

    mod = 91120
    MH   = 1
    no = 10
    if len(sys.argv) == 3:
        mod = int(sys.argv[1])
        MH  = int(sys.argv[2])

    # load initial values
    iA, idT, itau1, itau2, iC = loadInit(mod, MH)
    print("Init value tau1 = %.2e" %itau1)
    print("Init value tau2 = %.2e" %itau2)
    print("Init value A = %.2f" %(iA/iC))


    # declare variables...
    Tmin, Tmax = 0.003, 0.08     # Fitting time range
    time   = ROOT.RooRealVar("time", "time", Tmin, Tmax)
    time.setRange('nllRange', Tmin, Tmax)
    
    tau1 = ROOT.RooRealVar("tau1", "tau1", itau1, 0, 10)
    tau2 = ROOT.RooRealVar("tau2", "tau2", itau2, 0, 10)
    dT   = ROOT.RooRealVar("dT", "dT", 3.5e-3, -0.01, 0.01)
    A    = ROOT.RooRealVar("A", "A", iA, 0, 1000)


    # build Fitting P.D.F
    
    # build 1st exponential p.d.f
    model = ROOT.RooGenericPdf("model", "time spectra", "(1-exp(-(time)/tau1))*(A*exp(-(time)/tau2)+1)", ROOT.RooArgSet(time, tau1, A, tau2))


    # Load Data from toyMC
    filename = "/junofs/users/miaoyu/supernova/production/Data/5kpc/data_mod%d_cha1_mh%d_val.root"%(mod, MH)
    #filename = "/junofs/users/miaoyu/supernova/production/Data/2kpc/data_mod%d_cha1_mh%d_nuMass0.0_test%d.root"%(mod, MH, no)
    print(filename)
    dataFile  = ROOT.TFile(filename, 'read')
    inputTree = dataFile.Get('evt')

    evtID    = ROOT.RooRealVar("evtID", "evtID", 0, 500)

    newArgSet = ROOT.RooArgSet(time)

    #pdffilename = "test.pdf"
    filename = "/junofs/users/miaoyu/energy_model/fyyitter/energyModel_Fit/new_fitter/data/gamma/spectrum/"
    pdffilename = "./FitPanel_mod%d_mh%d_test%d_5kpc.pdf" %(mod, MH, no)
    can = ROOT.TCanvas("can", "can", 1200, 800)
    can.Print(pdffilename + "[")

    
    tau1_arr, tau2_arr, Nsig_arr, A_arr = [], [], [], []
    status_arr = []


    nStat = 500
    for iSub in range(1, nStat+1):
        print("Process toy dataset %d" %iSub)

        preSelData1D = ROOT.RooDataSet("dataGroup%d"%(iSub), "dataset with (e,t)", 
                       inputTree, ROOT.RooArgSet(time, evtID), 
                       "evtID==%d"%(iSub) )

    
        nFitEvts1D = round(preSelData1D.sumEntries())
        print('preSelData1D.sumEntries(): ', nFitEvts1D )

        fitdata = ROOT.RooDataSet('fitdata', 'fitdata', newArgSet)

        for iEvt in range(nFitEvts1D):
            argSet = preSelData1D.get(iEvt)
            timeTmp = argSet.getRealValue('time')
            weight = preSelData1D.weight()
            weightError = preSelData1D.weightError()
            
            newArgSet.setRealValue('time', timeTmp)

            fitdata.add( newArgSet, weight, weightError)

        #res = model.fitTo(fitdata, ROOT.RooFit.Range('nllRange'),
        #                ROOT.RooFit.Save(), ROOT.RooFit.Extended(),
        #                ROOT.RooFit.PrintLevel(-1))
        Nsig = ROOT.RooRealVar("Nsig", "Nsig", nFitEvts1D, 0, 1000)
        modelext = ROOT.RooExtendPdf("modelext", "modelext", model, Nsig)
        res = modelext.fitTo(fitdata, ROOT.RooFit.Range('nllRange'), 
                             ROOT.RooFit.Save(), ROOT.RooFit.Extended(),
                             ROOT.RooFit.PrintLevel(-1))
        nllVal = res.minNll()
        bestStatus = res.status()
        print("FIIIIIIIIT", tau1.getVal(), tau2.getVal(), A.getVal(), Nsig.getVal(), nllVal, bestStatus)

        if math.isnan(nllVal) != True :
            tau1_arr.append(tau1.getVal())
            tau2_arr.append(tau2.getVal())
            Nsig_arr.append(Nsig.getVal())
            A_arr.append(A.getVal())
            status_arr.append(bestStatus)

        can.cd()
        tframe = time.frame(ROOT.RooFit.Title("dataset fitting"), ROOT.RooFit.Range('nllRange'), ROOT.RooFit.Bins(30))
        fitdata.plotOn(tframe)
        modelext.plotOn(tframe, ROOT.RooFit.Range('nllRange'),
                        ROOT.RooFit.Normalization(nFitEvts1D, ROOT.RooAbsReal.NumEvent),
                        ROOT.RooFit.LineColor(2))
        tframe.Draw()
        can.Print(pdffilename + "[")
        del fitdata

    can.Print(pdffilename + "]")


    #with uproot3.recreate("test.root") as f:
    with uproot3.recreate("fitRes_mod%d_MH%d_test%d_5kpc.root"%(mod, MH, no)) as f:
        f["fit"] = uproot3.newtree({"tau1":"float64", "tau2":"float64", "A":"float64", "Nsig":"float64", "status":"int32"})
        f["fit"].extend({"tau1":tau1_arr, "tau2":tau2_arr, "A":A_arr, "Nsig":Nsig_arr, "status":status_arr})
    

