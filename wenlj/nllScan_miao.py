import numpy as np
import matplotlib.pyplot as plt
import ROOT


def getValue(strline):
    values = strline.split()
    x1 = float(values[1])
    x2 = float(values[2])
    x3 = float(values[3])
    return (x1, x2, x3)

def func(x, a0, t0, a1):
    return a0*np.power(x-t0, 2) + a1

def getNLL(filename):
    deltaT, nllVal, nsig = [], [], []
    file = open(filename, 'r')
    lines = file.readlines()
    count = 0
    
    tmp_deltaT, tmp_nllVal, tmp_nsig = [], [], []
    for k in range(len(lines)):
        if count == 25:
            deltaT.append(tmp_deltaT)
            nllVal.append(tmp_nllVal)
            nsig.append(tmp_nsig)
            count = 0
            tmp_deltaT, tmp_nllVal, tmp_nsig = [], [], []
        (dT, val, nEvt) = getValue( lines[k] )
        count += 1
        tmp_deltaT.append(dT)
        tmp_nllVal.append(val)
        tmp_nsig.append(nEvt)
        
        #nsig.append(nEvt)
    file.close()

    deltaT = np.array(deltaT)
    nllVal = np.array(nllVal)
    nsig = np.array(nsig)
    #locMinNLL = np.min(nllVal)

    return deltaT, nsig, nllVal


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







import ROOT
import namespace as ns
def loadPdf_Dataset(pdfname, pdfHistname, dataname, evtid, deltaT, pdffilename, can):
    
    # load PDF
    pdfFile = ROOT.TFile(pdfname, "read")
    inputHist = pdfFile.Get(pdfHistname)       # 2D input signal histogram

    Tmin = inputHist.GetXaxis().GetXmin()
    Tmax = inputHist.GetXaxis().GetXmax()
    Emin = inputHist.GetYaxis().GetXmin()
    Emax = inputHist.GetYaxis().GetXmax()
    print('fitting PDF, Tmin, Tmax, Emin, Emax: ', Tmin, Tmax, Emin, Emax)

    nuEnergy = ROOT.RooRealVar("nuEnergy", "nuEnergy", Emin, Emax)
    nuTime = ROOT.RooRealVar("nuTime", "nuTime", Tmin, Tmax)
    nuEnergy.setRange("fullRange", Emin, Emax)
    nuTime.setRange('fullRange', Tmin, Tmax)
    NEVENTS = getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")
    print('NEVENTS in full range: ', NEVENTS)

    binEthr = inputHist.GetYaxis().FindBin(0.1)
    print('Ethr, binEthr, lower bin edge of binEthr: ', 
          Ethr, ',', binEthr, ', ', inputHist.GetYaxis().GetBinLowEdge(binEthr)) 


    # 1D histogram in time domain
    ################# Input 1D Signal/Background 1D Histograms
    inputHist1D = inputHist.ProjectionX('_pdfsigpx', binEthr, -1, '')

    sigHist1D = ROOT.RooDataHist('sigHist1D', 'sigHist1D', 
          ROOT.RooArgList(nuTime), 
          ROOT.RooFit.Import(inputHist1D, True) )

    ## fitPdf1D from nllFit.py
    sigPdf1D = ROOT.RooHistPdf('sigPdf1D', 'sigPdf1D', ROOT.RooArgList(nuTime), sigHist1D, 2)
    nsig     = ROOT.RooRealVar("nsig","#signal events", NEVENTS, 0., NEVENTS*1.2)
    fitPdf1D = ROOT.RooExtendPdf("fitPdf1D", "fitPdf1D", sigPdf1D, nsig)

    # load dataset 

    dataFile = ROOT.TFile(dataname, "read") 
    inputTree = dataFile.Get('tFixedStat')

    Tmin, Tmax, Emin, Emax = -0.02, 0.025, 0.1, 4.0

    nuTime1D = ROOT.RooRealVar("nuTime1D", "nuTime1D", Tmin, Tmax)
    evtID = ROOT.RooRealVar("evtID", "evtID", 0, 5000)

    newArgSet = ROOT.RooArgSet(nuTime)

    preSelData1D = ROOT.RooDataSet("dataGroup%d"%(evtid), "dataset with (e,t)", 
                    inputTree, ROOT.RooArgSet(nuTime1D, evtID), 
                    'evtID==%d' %(evtid))
    nFitEvts1D = round(preSelData1D.sumEntries())


    fitdata = ROOT.RooDataSet('fitdata', 'fitdata', newArgSet)
    nuTime.setRange('nllRange', Tmin + deltaT, Tmax + deltaT)
    for iEvt in range(nFitEvts1D):
        argSet = preSelData1D.get(iEvt)
        timeTmp = argSet.getRealValue('nuTime1D')
        weight = preSelData1D.weight()
        weightError = preSelData1D.weightError()
        newArgSet.setRealValue('nuTime', timeTmp + deltaT)
        fitdata.add( newArgSet, weight, weightError)

    ## Drawing

    xframe = nuTime.frame(ROOT.RooFit.Title("data group %d, dT: %8.5f"%(evtid, deltaT)), 
                    ROOT.RooFit.Range('nllRange'), ROOT.RooFit.Bins(100))
    xframe.GetYaxis().SetTitleOffset(1.4)
    can.cd()
    fitdata.plotOn(xframe)
    fitPdf1D.plotOn(xframe, ROOT.RooFit.Range('nllRange'),
            ROOT.RooFit.Normalization(nFitEvts1D, ROOT.RooAbsReal.NumEvent),
            ROOT.RooFit.LineColor(2))
    xframe.Draw()

    can.Print(pdffilename) 


def drawRes(can, mass, nll, deltaT, nsig, pdffilename) :

    gr0 = ROOT.TGraph()
    gr1 = ROOT.TGraph()
    gr2 = ROOT.TGraph()

    for i in range(20):
        gr0.SetPoint(i, mass[i], nll[i])
        gr1.SetPoint(i, mass[i], deltaT[i])
        gr2.SetPoint(i, mass[i], nsig[i])


    gr0.SetLineColor(9)
    gr0.SetLineWidth(2)
    gr0.SetMarkerColor(9)
    gr0.SetMarkerStyle(20)
    gr1.SetLineColor(9)
    gr1.SetLineWidth(2)
    gr1.SetMarkerColor(9)
    gr1.SetMarkerStyle(20)
    gr2.SetLineColor(9)
    gr2.SetLineWidth(2)
    gr2.SetMarkerColor(9)
    gr2.SetMarkerStyle(20)

    can.cd(1)
    gr0.Draw("APC")
    can.cd(2)
    gr1.Draw("APC")
    can.cd(3)
    gr2.Draw("APC")

    can.Print(pdffilename)



if __name__ == "__main__" :

    nuMass1 = 0.0
    #nuMass2 = 0.2
    modelNum = 82503
    chaname = 1
    chaName = ['nup','nue','IBD']
    dist = 10.0
    MH = 2
    Ethr = 0.1
    group = 1
    mNu, minNLL = [], []
    prefix = "./dataset/%2.1fkpc/TEvisDATA_" %(dist)


    minVal_nuMass = np.zeros((500, 20))
    deltaT_nuMass = np.zeros((500, 20))
    nsig_nuMass = np.zeros((500, 20))
    haveBadMass = np.zeros(500)
    
    for k in range(0, 20):
        nuMass2 = k * 0.1
        fname = prefix+'mod%d_cha%d%s_dataMH%d_fitMH%d_data%.1feV_pdf%.1feV_%2.1fkpc_Ethr0.1MeV_group1_num2_1DFitRaw.txt'%(modelNum, chaname, chaName[chaname], MH, MH, nuMass1, nuMass2, dist)
        deltaT, nsig, nllVal = getNLL(fname)

        print(fname)

        # Drop all events with abnormal fitting points
        
        evtid = 0
        for tarr, narr, nllarr in zip(deltaT, nsig, nllVal):
            locMinNll, locMaxNll = [], []
            for y in range(1, len(nllarr)-1):
                if nllarr[y-1] < nllarr[y] and nllarr[y+1] < nllarr[y] :
                    locMaxNll.append(nllarr[y])
                if nllarr[y-1] > nllarr[y] and nllarr[y+1] > nllarr[y] :
                    locMinNll.append(nllarr[y])
            
            if len(locMinNll) > 1 or len(locMaxNll) > 0 :
                #print(tarr)
                #print(nllarr)
                #mask.append(evtid)
                
                #minVal_nuMass[evtid, k] = 999
                #deltaT_nuMass[evtid, k] = 999
                #nsig_nuMass[evtid, k]   = 999
                minVal_nuMass[evtid, k] = np.min(nllarr)
                deltaT_nuMass[evtid, k] = tarr[nllarr.argmin()]
                nsig_nuMass[evtid, k]   = narr[nllarr.argmin()]
                haveBadMass[evtid]      = 1


            else :    # normal case
                minVal_nuMass[evtid, k] = np.min(nllarr)
                deltaT_nuMass[evtid, k] = tarr[nllarr.argmin()]
                nsig_nuMass[evtid, k]   = narr[nllarr.argmin()]

            evtid += 1

            if evtid >= 500:
                break


    # Plotting ...
    goodPlot = 0
    goodSelf = 0
    plotNum = 0
    nuMass = np.arange(0, 2, 0.1)

    #pdfilename = ns.fitResPdfName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)
    can = ROOT.TCanvas("can", "can", 1200, 800)
    pdffilename = "./FittingVisal_nllDeltaTBad_mh2.pdf"
    #can.Print(pdffilename + '[')

    can2 = ROOT.TCanvas("can2", "can2", 1200, 800)
    can2.Divide(3, 1)
    can2.cd(1).SetLeftMargin(0.15)
    can2.cd(1).SetRightMargin(0.05)
    can2.cd(2).SetLeftMargin(0.15)
    can2.cd(2).SetRightMargin(0.05)
    can2.cd(3).SetRightMargin(0.15)

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(12, 3))

    for i in range(500):

        ### check if very sharp changes exist ...
        #if minVal_nuMass[i, 4] - minVal_nuMass[i, 0] > 2:
        #    print("Sharp EventId : %d" %i)

        if minVal_nuMass[i, 0] - np.min(minVal_nuMass[i, :]) > 100 :
            continue

        ax0.plot(nuMass, minVal_nuMass[i, :]-np.min(minVal_nuMass[i, :]), "-")
        ax1.plot(nuMass, deltaT_nuMass[i, :], "-")
        ax2.plot(nuMass, nsig_nuMass[i, :] - nsig_nuMass[i, 0], "-")

        """

        if haveBadMass[i] == 1:    # some nuMass fitting results in this event, nll-deltaT is not a good shape
            if plotNum < 10:
            #    ax0.plot(nuMass, minVal_nuMass[i, :]-np.min(minVal_nuMass[i, :]), "-")
            #    ax1.plot(nuMass, deltaT_nuMass[i, :], "-")
            #    ax2.plot(nuMass, nsig_nuMass[i, :] - nsig_nuMass[i, 0], "-")

    
                # Draw Bad Fitting details :
                for mm in np.arange(0, 2.0, 0.1):
                    print("*************************************** plotNum %d -> nuMass %.1f **************************************" %(plotNum, mm))
                    pdfname = ns.pdfOldFileName(modelNum, dist, chaname, mm, MH)
                    pdfHistname = ns.pdfEvisHist2DName(modelNum, chaname, MH)
                    dataname = ns.dataFileName(modelNum, dist, chaname, 0, MH, 0.1, 1)

                    loadPdf_Dataset(pdfname, pdfHistname, dataname, i, deltaT_nuMass[i, int(mm*10)], pdffilename, can)

                drawRes(can2, nuMass, minVal_nuMass[i, :], deltaT_nuMass[i, :], nsig_nuMass[i, :], pdffilename)

                plotNum += 1
            continue

        else :

            goodSelf += 1
            # One more shape checks on NLL profile :

            locMinNll, locMaxNll = 0, 0
            for j in range(1, 19, 1):
                if minVal_nuMass[i, j] > minVal_nuMass[i, j-1] and minVal_nuMass[i, j] > minVal_nuMass[i, j+1]:
                    locMaxNll += 1
                if minVal_nuMass[i, j] < minVal_nuMass[i, j-1] and minVal_nuMass[i, j] < minVal_nuMass[i, j+1]:
                    locMinNll += 1

            if locMinNll > 1 or locMaxNll > 0 :
                # Draw bad fitting cases
                #if plotNum < 10:
                ##    ax0.plot(nuMass, minVal_nuMass[i, :]-np.min(minVal_nuMass[i, :]), "-")
                ##    ax1.plot(nuMass, deltaT_nuMass[i, :], "-")
                ##    ax2.plot(nuMass, nsig_nuMass[i, :] - nsig_nuMass[i, 0], "-")

    
                #    # Draw Bad Fitting details :
                #    for mm in np.arange(0, 2.0, 0.1):
                #        print("*************************************** plotNum %d -> nuMass %.1f **************************************" %(plotNum, mm))
                #        pdfname = ns.pdfOldFileName(modelNum, dist, chaname, mm, MH)
                #        pdfHistname = ns.pdfEvisHist2DName(modelNum, chaname, MH)
                #        dataname = ns.dataFileName(modelNum, dist, chaname, 0, MH, 0.1, 1)

                #        loadPdf_Dataset(pdfname, pdfHistname, dataname, i, deltaT_nuMass[i, int(mm*10)], pdffilename, can)

                #    drawRes(can2, nuMass, minVal_nuMass[i, :], deltaT_nuMass[i, :], nsig_nuMass[i, :], pdffilename)

                #    plotNum += 1


                continue
    

            #plt.plot(nuMass, minVal_nuMass[i, :]-np.min(minVal_nuMass[i, :]), "-")

            #if plotNum < 4:
            #    ax0.plot(nuMass, minVal_nuMass[i, :]-np.min(minVal_nuMass[i, :]), "-")
            #    ax1.plot(nuMass, deltaT_nuMass[i, :], "-")
            #    ax2.plot(nuMass, nsig_nuMass[i, :] - nsig_nuMass[i, 0], "-")
            #    plotNum += 1

            goodPlot += 1
        """

    #print("good fitting for all nuMass in one event number : %d " %goodSelf)
    #print("Final good plots number = %d" %(goodPlot) )
    #plt.xlabel("fitNuMass/eV")
    #plt.ylabel("nll")
    
    ax0.set_xlabel("fitNuMass/eV")
    ax0.set_ylabel("nll")
    #
    ax1.set_xlabel("fitNuMass/eV")
    ax1.set_ylabel("deltaT/ns")

    ax2.set_xlabel("fitNuMass/eV")
    ax2.set_ylabel("nEvt")
    
    plt.show()


    #can.Print(pdffilename + ']')



