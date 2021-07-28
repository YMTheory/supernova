import ROOT
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import json
import scipy
import numpy as np
import gc
import time
from array import array
import namespace as ns

def setGraphStyle(gr, markerStyle, lineStyle):
    gr.SetMarkerStyle(markerStyle[0])
    gr.SetMarkerSize(markerStyle[1])
    gr.SetMarkerColor(markerStyle[2])
    gr.SetLineStyle(lineStyle[0])
    gr.SetLineWidth(lineStyle[1])
    gr.SetLineColor(lineStyle[2])

def getIntegral(hist2d, x1, x2, y1, y2, opt):
    binx1 = hist2d.GetXaxis().FindBin(x1)
    binx2 = hist2d.GetXaxis().FindBin(x2)
    biny1 = hist2d.GetYaxis().FindBin(y1)
    biny2 = hist2d.GetYaxis().FindBin(y2)

    return hist2d.Integral(binx1, binx2, biny1, biny2, opt)


def getValue(strline):
    values = strline.split()
    x1 = float(values[0])
    x2 = float(values[1])
    x3 = float(values[2])
    return (x1, x2, x3)

def getNLL(filename):
    deltaT, nllVal, nsig = [], [], []
    file = open(filename, 'r')
    lines = file.readlines()
    for k in range(len(lines)):
        (dT, val, nEvt) = getValue( lines[k] )
        deltaT.append(dT)
        nllVal.append(val)
        nsig.append(nEvt)
    file.close()

    deltaT = np.array(deltaT)
    nllVal = np.array(nllVal)
    print(deltaT)
    print(nllVal)
    locMinNLL = np.min(nllVal)

    return deltaT, nllVal, locMinNLL

if __name__ == "__main__":

    if len(sys.argv) < 11:
        print("Usage: %s [modelNum] [channel] [nuType] [nuMass1] [dataMH] [nuMass2] [pdfMH] [dist] [Ethr] [group]"%(sys.argv[0]))
        print("modelNum: check the SNsim/simulation/data directory")
        print("channel : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("nuMass1 : neutrino mass that used to generate dataset, in eV unit")
        print("dataMH  : NMO that used to generate dataset, 0 for NoOsc, 1 for NH, 2 for IH")
        print("nuMass2 : neutrino mass that used in PDF to fit the dataset, in eV unit")
        print("pdfMH   : NMO that used to generate PDF")
        print("dist    : SN distance in kpc unit")
        print("Ethr    : fit Energy threshold")
        print("group   : index of the groups")
        print("example : python nllFit.py 82503 1 0 0.0 2 0.2 1 10.0 0.1 1")
        sys.exit(1)

    ##################################
    # Input variables
    ROOT.RooFit.PrintLevel(-1)
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setSilentMode(True)
    # choose SN model
    #modelNum = 82503
    modelNum = int( sys.argv[1] )

    chaName = ['nup','nue','IBD']
    #chaname = snDet.IBD
    chaname = int( sys.argv[2] )
    print("channelName: ", chaname, chaName[chaname])

    nuType = int( sys.argv[3] )
    # -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x
    nuTypeName = ['nu_e', 'anti_nue', 'nu_x']
    print("neutrino type: ", nuType, nuTypeName[nuType])
    
    nuMass1 = float( sys.argv[4] )# unit: eV
    print('nuMass for producing dataset : ', nuMass1, ' eV')

    MHName = ['NoOsc', 'NH', 'IH']
    dataMH = int( sys.argv[5] )
    print('NMO for producing dataset : ', dataMH, MHName[dataMH])

    nuMass2 = float( sys.argv[6] )# unit: eV
    print('nuMass for producing PDF: ', nuMass2, ' eV')

    fitMH = int( sys.argv[7] )
    print('NMO for producing PDF: ', fitMH, MHName[fitMH])

    dist = float( sys.argv[8] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[9] )# unit: kpc
    print('fit Ethr: ', Ethr, ' MeV')

    group = int( sys.argv[10] )# group number
    print('group number: ', group)




    ################################################################################33
    # # Obtain the fitting PDF
    filename = ns.pdfFileName(modelNum, dist, chaname, nuMass2, fitMH)
    print(filename)
    pdfFile  = ROOT.TFile(filename, 'read')
    inputHist = pdfFile.Get(ns.pdfSigHist2DName(modelNum, chaname, fitMH))
    inputBKG  = pdfFile.Get(ns.pdfBkgHist2DName(modelNum, chaname, fitMH))

    Tmin = inputHist.GetXaxis().GetXmin()
    Tmax = inputHist.GetXaxis().GetXmax()
    Emin = inputHist.GetYaxis().GetXmin()
    Emax = inputHist.GetYaxis().GetXmax()
    print('fitting PDF, Tmin, Tmax, Emin, Emax: ', Tmin, Tmax, Emin, Emax)
    nuEnergy = ROOT.RooRealVar("nuEnergy", "nuEnergy", Emin, Emax)
    nuTime   = ROOT.RooRealVar("nuTime", "nuTime", Tmin, Tmax)
    nuEnergy.setRange('fullRange', Emin, Emax)
    nuTime.setRange('fullRange', Tmin, Tmax)

    offset = ROOT.RooRealVar("offset", "offset", -0.01, 0.023)
    shiftNuT = ROOT.RooFormulaVar("shiftNuT", "shiftNuT", "@0 + @1",\
                            ROOT.RooArgList(nuTime, offset))

    # Integrated Number of Events in the full range of the histograms
    NEVENTS = getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")
    print('NEVENTS in full range: ', NEVENTS)

    #NEVENTS = NEVENTS * 4 # special test for 10 kpc, group # == 999
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    #NEVENTS = NEVENTS * 25 # special test for 10 kpc, group # == 88X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    #NEVENTS = NEVENTS * 100 # special test for 10 kpc, group # == 77X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    #NEVENTS = int(NEVENTS * 6.25) # special test for 5 kpc, group # == 88X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    # change the range for pre-selection to the dataset
    if chaname == 2:  # snDet.IBD
        Tmin, Tmax, Emin, Emax = -0.02, 0.2, 1.0, 10.0
    if chaname == 1:  # snDet.NuE
        Tmin, Tmax, Emin, Emax = -0.02, 0.025, 0.1, 10.0

    #if fitMH==0:
    #    Tmin, Tmax = -0.015, 0.02
    #if fitMH==1:
    #    Tmin, Tmax = -0.03, 0.02
    #    #Tmin, Tmax = -0.04, 0.08
    #if fitMH==2:
    #    Tmin, Tmax = -0.03, 0.02
        #Tmin, Tmax = -0.04, 0.08

    nuTime.setRange('nllRange', Tmin, Tmax)

    # 1D histogram in time domain
    binEthr = inputHist.GetYaxis().FindBin(Ethr)
    print('Ethr, binEthr, lower bin edge of binEthr: ', 
          Ethr, ',', binEthr, ', ', inputHist.GetYaxis().GetBinLowEdge(binEthr)) 

    binEvismax = inputHist.GetYaxis().FindBin(Emax)
    print('Ethr, binEvismax, lower bin edge of binEvismax: ', 
          Ethr, ',', binEvismax, ', ', inputHist.GetYaxis().GetBinLowEdge(binEvismax)) 

    inputHist1D = inputHist.ProjectionX('_pdfsigpx', binEthr, binEvismax, '')
    inputBKG1D  =  inputBKG.ProjectionX('_pdfbkgpx', binEthr, binEvismax, '')

    sigHist = ROOT.RooDataHist('sigHist', 'sigHist', 
          ROOT.RooArgList(nuTime, nuEnergy), 
          ROOT.RooFit.Import(inputHist, True) )

    #bkgHist = ROOT.RooDataHist('bkgHist', 'bkgHist', 
    #      ROOT.RooArgList(nuTime, nuEnergy), 
    #      ROOT.RooFit.Import(inputBKG, True) )

    sigPdf = ROOT.RooHistPdf('sigPdf', 'sigPdf', 
          ROOT.RooArgList(nuTime, nuEnergy), sigHist, 2)

    # 1D PDF in time domain
    sigHist1D = ROOT.RooDataHist('sigHist1D', 'sigHist1D', 
          ROOT.RooArgList(nuTime), 
          ROOT.RooFit.Import(inputHist1D, True) )

    sigPdf1D = ROOT.RooHistPdf('sigPdf1D', 'sigPdf1D', 
           ROOT.RooArgList(shiftNuT), ROOT.RooArgList(nuTime), sigHist1D, 1)

    nsig = ROOT.RooRealVar("nsig","#signal events", NEVENTS, 0., NEVENTS*10.0)
    fitPdf1D = ROOT.RooExtendPdf('fitPdf1D', 'fitPdf1D', sigPdf1D, nsig)

    print('pdf integral: ', \
         fitPdf1D.analyticalIntegralWN(0, ROOT.RooArgSet(nuTime), 'nllRange'))
    print('full integral: ', \
         fitPdf1D.analyticalIntegralWN(0, ROOT.RooArgSet(nuTime), 'fullRange'))

    ##################################
    # Obtain data sets from the root files
    filename = ns.dataFileName(modelNum, dist, chaname, nuMass1, dataMH, Ethr, group)
    print(filename)
    dataFile  = ROOT.TFile(filename, 'read')
    inputTree = dataFile.Get('tFixedStat')

    nuTime1D = ROOT.RooRealVar("nuTime1D", "nuTime1D", Tmin, Tmax)
    evtID    = ROOT.RooRealVar("evtID", "evtID", 0, 1000)

    resFiles = ns.fitResFileName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)
    #resfileRaw = open(resFiles[0], "w")
    resfileSummary = open(resFiles[1], 'w')

    can = ROOT.TCanvas("can", "can", 1200, 800)
    can.SetTopMargin(0.1)
    can.SetLeftMargin(0.1)
    ##pdffilename = prefix+'_data%2.1feV_mh%d_pdf%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_1DFit_Tmax%3.2f.pdf'%(nuMass1, fitMH, nuMass2, dist, Ethr, group, Tmax)
    pdffilename = ns.fitResPdfName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)
    can.Print(pdffilename + '[')

    can2 = ROOT.TCanvas("can2", "can2", 1200, 1000)
    can2.Divide(2,2)
    can2.cd(1).SetLeftMargin(0.15)
    can2.cd(1).SetRightMargin(0.05)
    can2.cd(2).SetRightMargin(0.05)

    nStat = 500
    for iSub in range(1,nStat+1):

        print('Process toy dataset %d'%iSub)
        preSelData1D = ROOT.RooDataSet("dataGroup%d"%(iSub), "dataset with (e,t)", 
                        inputTree, ROOT.RooArgSet(nuTime1D, evtID), 
                        "evtID==%d"%(iSub) )
    
        nFitEvts1D = round(preSelData1D.sumEntries())
        print('preSelData1D.sumEntries(): ', preSelData1D.sumEntries())
        print('preSelData1D.sumEntries(): ', nFitEvts1D )

        ################# Likelihood fit
        nllList, offsetList, nsigList = [], [], []
        #dTexp =  nuMass2*0.005
        #offset.setVal(dTexp)

        ################## New Data Set
        newArgSet = ROOT.RooArgSet(nuTime)
        fitdata = ROOT.RooDataSet('fitdata', 'fitdata', newArgSet)
        
        for iEvt in range(nFitEvts1D):
            argSet = preSelData1D.get(iEvt)
            timeTmp = argSet.getRealValue('nuTime1D')
            #energyTmp = argSet.getRealValue('nuEnergy')
            if timeTmp>Tmax or timeTmp<Tmin:
                print('timeTmp>Tmax or timeTmp<Tmin')
                continue

            weight = preSelData1D.weight()
            weightError = preSelData1D.weightError()
            if iEvt==0:
                print('time, weight, weightError: ', 
                        timeTmp, weight, weightError)

            newArgSet.setRealValue('nuTime', timeTmp)
            #newArgSet.setRealValue('nuEnergy', energyTmp)
            fitdata.add( newArgSet, weight, weightError)

        ################## Fit
        nFitEvts1D = round(fitdata.sumEntries())
        print('init nsig: ', nsig.getVal())

        offsetDist = 0.

        if dataMH == 2:
            if dist==10.0:
                if nuMass2<1.5:
                    offsetDist = 0.003 * (nuMass2 - nuMass1)
                else:
                    offsetDist = 0.005 * (nuMass2 - nuMass1)

                if fitMH == 1:
                    offsetDist = 0.01 + offsetDist

            if dist==5.0:
                offsetDist = 0.0025 * (nuMass2 - nuMass1)
                if fitMH == 1:
                    offsetDist = 0.01 + offsetDist

            if dist==3.0 or dist==2.0:
                offsetDist = 0.0015 * (nuMass2 - nuMass1)
                if fitMH == 1:
                    offsetDist = 0.012 + offsetDist

        if dataMH == 1:
            if dist==10.0:
                if nuMass2<1.5:
                    offsetDist = 0.003 * (nuMass2 - nuMass1)
                else:
                    offsetDist = 0.005 * (nuMass2 - nuMass1)

                if fitMH == 2:
                    offsetDist = offsetDist - 0.005

            if dist==5.0:
                offsetDist = 0.0025 * (nuMass2 - nuMass1)
                if fitMH == 2:
                    offsetDist = offsetDist - 0.005

            if dist==3.0 or dist==2.0:
                offsetDist = 0.0018 * (nuMass2 - nuMass1)
                if fitMH == 2:
                    offsetDist = offsetDist - 0.005

        bestFit = []
        tScan = []
        nTry = 10
        for iFit in range(nTry):
            dTexp = (2*iFit+1-nTry)*0.001/4 + offsetDist
            tScan.append(dTexp)

            offset.setVal(dTexp)

            res = fitPdf1D.fitTo(fitdata, ROOT.RooFit.Range('nllRange'),
                            ROOT.RooFit.Save(), #ROOT.RooFit.Extended(),
                            ROOT.RooFit.Minimizer('Minuit2'),
                            ROOT.RooFit.PrintLevel(-1))
            if res.status()!=0:
                res = fitPdf1D.fitTo(fitdata, ROOT.RooFit.Range('nllRange'),
                                ROOT.RooFit.Save(), #ROOT.RooFit.Extended(),
                                ROOT.RooFit.Minimizer('Minuit2'),
                                ROOT.RooFit.PrintLevel(-1))
                
            bestFit.append( [res.minNll(), nsig.getVal(), \
                             offset.getVal(), res.status()] )

            nllList.append(res.minNll())
            offsetList.append(offset.getVal())
            nsigList.append(nsig.getVal())
        print('Good fit num: ', len(bestFit))
        #print(bestFit, sep=', ', end='\n')
        resArr = sorted(bestFit, key=lambda x:x[0], reverse=False)
        print(resArr, sep=', ', end='\n')
        fitRes = resArr[0]
        print('nllVal, nsig, offset: %10.3f, %10.3f, %10.6f, status: %d'%\
              (fitRes[0], fitRes[1], fitRes[2], fitRes[3]))
        #      (nllVal, nsig.getVal(), offset.getVal(), res.status()))
        #nsig.setVal(295.868)
        #offset.setVal(0.006)
        #print('nsig, offset: %10.3f, %10.6f'%\
        #      (nsig.getVal(), offset.getVal()))

        if iSub % 10 == 0:

            xframe = nuTime.frame(ROOT.RooFit.Title("dataGroup%d"%iSub), 
                            ROOT.RooFit.Range('nllRange'), ROOT.RooFit.Bins(100))
            xframe.GetYaxis().SetTitleOffset(0.8)
            xframe.GetXaxis().SetTitle('nuTime (s) at dataset %d'%iSub)
            can.cd()
            dataPlot = fitdata.plotOn(xframe)
            # fitPdf.plotOn(xframe, ROOT.RooFit.ProjectionRange('nllRange'),
            #                 ROOT.RooFit.LineColor(2) )
                                #ROOT.RooFit.DrawOption("E"),
                                #ROOT.RooFit.FillColor(ROOT.kOrange),
                                #ROOT.RooFit.MoveToBack() )
            pdfPlot = fitPdf1D.plotOn(xframe, ROOT.RooFit.Range('nllRange'),
                    ROOT.RooFit.Normalization(nsig.getVal(), ROOT.RooAbsReal.NumEvent),
                    ROOT.RooFit.LineColor(2))
            xframe.Draw()
            can.Print(pdffilename)

            ###### Draw the nll profile of iSub toy dataset -> 
            tScanX    = array('d', tScan)
            nllListY  = array('d', nllList)
            nsigListY = array('d', nsigList) 
            offsetListX = array('d', offsetList)
            gr1 = ROOT.TGraph(len(nllList), tScanX, nllListY)
            gr2 = ROOT.TGraph(len(nllList), tScanX, nsigListY)
            gr3 = ROOT.TGraph(len(nllList), offsetListX, nsigListY)
            setGraphStyle(gr1, [20, 1.2, 1], [2, 1, 4])
            setGraphStyle(gr2, [20, 1.2, 1], [2, 1, 4])
            setGraphStyle(gr3, [20, 1.2, 1], [2, 1, 4])
            can2.cd(1)
            gr1.Draw('APC')
            gr1.GetXaxis().SetTitle('#delta_{T} (s)')
            gr1.GetYaxis().SetTitle('nll')
            can2.cd(2)
            gr2.Draw('APC')
            gr2.GetXaxis().SetTitle('#delta_{T} (s)')
            gr2.GetYaxis().SetTitle('number of events')
            can2.cd(3)
            gr3.Draw('APC')
            gr3.GetXaxis().SetTitle('#offset (s)')
            gr3.GetYaxis().SetTitle('nll')
            can2.Print(pdffilename)

        fitdata.Clear()

        # iSub, offset, nll, nsig, nEvtTrue, status
        #for k in range(len(nllList)):
        #    resfileRaw.write('%d  %9.6f  %12.3f  %12.3f\n'%(iSub,offsetList[k], nllList[k], nsigList[k]))
        resfileSummary.write('%d  %9.6f  %12.3f  %12.3f  %d  %d\n'%(iSub, fitRes[2], fitRes[0], fitRes[1], nFitEvts1D, fitRes[3]))

    resfileSummary.close()
    pdfFile.Close()
    can.Print(pdffilename+']')
