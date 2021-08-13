#############################################################
### In this script, I try to do 2D(time, energy) fitting  ### 
### based on RooFit.                                      ###
### DeltaT scanning runs twice, one coarse and one fine   ###
### with the initial values from first run.               ###
### Outputs for all fitting results.                      ###
###                                                       ### 
###                 miaoyu@ihep.ac.cn                     ### 
#############################################################


import ROOT
import math
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import sys
import json
import scipy
import numpy as np
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

    # print('LowEdge(binx1)',  hist2d.GetXaxis().GetBinLowEdge(binx1))
    # print('LowEdge(binx2+1)',  hist2d.GetXaxis().GetBinLowEdge(binx2+1))
    # print('LowEdge(biny1)',  hist2d.GetYaxis().GetBinLowEdge(biny1))
    # print('LowEdge(biny2+1)',  hist2d.GetYaxis().GetBinLowEdge(biny2+1))

    return hist2d.Integral(binx1, binx2, biny1, biny2, opt)

# def plotDataSet(xList, yList):
#     plotx = np.arange(-0.2, 0.2)
#     plt.plot(xList, yList,  "o", ms=5, label="Sinnock 1969")
#     plt.xlabel('$\delta_{T}$ (s)')
#     plt.ylabel('NLL')
#     titlename = 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)
#     titlename += '_data%2.1feV_pdf%2.1feV'%(nuMass1, nuMass2)
#     plt.title(titlename)
#     pdf.savefig()

# def printHist(inputHist):
#     Tmin, Tmax, Emin, Emax = -0.05, 0.3, 0.2, 5.0
#     print('Integral (%4.2f, %4.2f), (%4.2f, %4.2f): %10.3f'
#         %(Tmin, Tmax, Emin, Emax, getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")*NEVENTS))

#     Tmin, Tmax, Emin, Emax = -0.05, 0.2, 0.2, 3.0
#     print('Integral (%4.2f, %4.2f), (%4.2f, %4.2f): %10.3f'
#         %(Tmin, Tmax, Emin, Emax, getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")*NEVENTS))
    
#     Tmin, Tmax, Emin, Emax = -0.05, 0.3, 0.2, 3.0
#     print('Integral (%4.2f, %4.2f), (%4.2f, %4.2f): %10.3f'
#         %(Tmin, Tmax, Emin, Emax, getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")*NEVENTS))

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


    # ###################################
    # # Obtain the fitting PDF
    #prefix = "./etSpec/fineSpec/TEvisPDF_"
    #prefix += 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],fitMH)
    #filename = prefix+'_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25.root'%(nuMass2, dist)
    filename = ns.pdfOldFileName(modelNum, dist, chaname, nuMass2, fitMH)
    print(filename)
    pdfFile  = ROOT.TFile(filename, 'read')
    inputHist = pdfFile.Get(ns.pdfEvisHist2DName(modelNum, chaname, fitMH))
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

    # Integrated Number of Events in the full range of the histograms
    NEVENTS = getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")
    print('NEVENTS in full range: ', NEVENTS)

    #NEVENTS = NEVENTS * 4 # special test for 10 kpc, group # == 999
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    #NEVENTS = NEVENTS * 25 # special test for 10 kpc, group # == 88X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    NEVENTS = NEVENTS * 100 # special test for 10 kpc, group # == 77X
    print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    # 1D histogram in time domain
    binEthr = inputHist.GetYaxis().FindBin(Ethr)
    print('Ethr, binEthr, lower bin edge of binEthr: ', 
          Ethr, ',', binEthr, ', ', inputHist.GetYaxis().GetBinLowEdge(binEthr)) 

    ################# Input 1D Signal/Background 1D Histograms

    sigHist = ROOT.RooDataHist('sigHist', 'sigHist', 
          ROOT.RooArgList(nuTime, nuEnergy), 
          ROOT.RooFit.Import(inputHist, True) )

    bkgHist = ROOT.RooDataHist('bkgHist', 'bkgHist', 
          ROOT.RooArgList(nuTime, nuEnergy), 
          ROOT.RooFit.Import(inputBKG, True) )

    sigPdf = ROOT.RooHistPdf('sigPdf', 'sigPdf', 
          ROOT.RooArgList(nuTime, nuEnergy), sigHist, 2)

    bkgPdf = ROOT.RooHistPdf('bkgPdf', 'bkgPdf', 
          ROOT.RooArgList(nuTime, nuEnergy), bkgHist, 2)


    #scaleFac = ROOT.RooRealVar('scaleFac', 'scaleFac', )
    nsig = ROOT.RooRealVar("nsig","#signal events", NEVENTS, 0., NEVENTS*1.2)
    nbkg = ROOT.RooRealVar("nbkg","#background events", 0.0)


    # fitPdf1D from nllFit.py
    fitPdf2D = ROOT.RooExtendPdf("fitPdf2D", "fitPdf2D", sigPdf, nsig)


    ##################################
    # Obtain data sets from the root files
    filename = ns.dataFileName(modelNum, dist, chaname, nuMass1, dataMH, Ethr, group)
    print(filename)
    dataFile  = ROOT.TFile(filename, 'read')
    inputTree = dataFile.Get('tFixedStat')

    # change the range for pre-selection to the dataset
    if chaname == 2:  # snDet.IBD
        Tmin, Tmax, Emin, Emax = -0.02, 0.2, 1.0, 10.0
    if chaname == 1:  # snDet.NuE
        Tmin, Tmax, Emin, Emax = -0.02, 0.025, 0.1, 10.0

    nuTime2D = ROOT.RooRealVar("nuTime2D", "nuTime2D", Tmin, Tmax)
    nuEnergy2D = ROOT.RooRealVar("nuEnergy2D", "nuEnergy2D", Emin, Emax)
    evtID    = ROOT.RooRealVar("evtID", "evtID", 0, 5000)


    resFiles = ns.fit2DResFileName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)
    resfileRaw = open(resFiles[0], "w")
    resfileSummary = open(resFiles[1], 'w')
    #print(resFiles[0])
    #print(resFiles[1])
    resfileRaw = open(resFiles[0], 'w')
    
    #can = ROOT.TCanvas("can", "can", 1200, 800)
    ##pdfilename = prefix+'_data%2.1feV_pdf%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_1DFit.pdf'%(nuMass1, nuMass2, dist, Ethr, group)
    #pdfilename = ns.fitResPdfName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)
    #can.Print(pdfilename + '[')

    #can2 = ROOT.TCanvas("can2", "can2", 1200, 800)
    #can2.Divide(2,1)
    #can2.cd(1).SetLeftMargin(0.15)
    #can2.cd(1).SetRightMargin(0.05)
    #can2.cd(2).SetRightMargin(0.05)

    newArgSet2D = ROOT.RooArgSet(nuTime2D, nuEnergy2D)

    nStat = 500
    for iSub in range(1, nStat+1):
        print('Process toy dataset %d'%iSub)
        # Get RooDataSet from the tree

        preSelData2D = ROOT.RooDataSet("dataGroup%d"%(iSub), "2d dataset with (e, t)", 
                        inputTree, ROOT.RooArgSet(nuTime2D, nuEnergy, evtID),
                        "evtID==%d"%(iSub))
    
        nFitEvts2D = round(preSelData2D.sumEntries())
        print('preSelData2D.sumEntries(): ', nFitEvts2D )

        ################# Likelihood Coarse scan
        #print('################ Start of Coase Scan ################')
        nllList, tScan, nsigList = [], [], []
        dTexp = 5.14e-3 * (nuMass2*nuMass2-nuMass1*nuMass1) * (100.0/9/9) * (dist/10.)
        #dTexp = 5.14e-3 * (nuMass2*nuMass2) * (100.0/10/10) * (dist/10.)
        #print('dTexp: ', dTexp)
        dTstep = 0.0005
        NScans = 35
        locMinNLL, bestfitDT = 100.0, 1.0
        extDTmax, extDTmin = 0.05, -0.05
        Tmin, Tmax = -0.015, 0.02

        for iScan in range(NScans):
            dT = dTexp + (2*iScan-NScans)*dTstep/2
            if dT>extDTmax or dT<extDTmin:
                print('Warning!: the data may be shifted outside the fitting range')
                continue

            ################## New Data Set
            #print('################ iSub %d, Start of Coase Scan %d ################'%(iSub, iScan))
            fitdata2D = ROOT.RooDataSet('fitdata2D', 'fitdata2D', newArgSet2D)
            # change the range for NLL fitting the dataset
            nuTime2D.setRange('nllRange', Tmin + dT, Tmax + dT)
            nuEnergy2D.setRange('nllRange', Emin, Emax)

            for iEvt in range(nFitEvts2D):
                argSet = preSelData2D.get(iEvt)
                timeTmp = argSet.getRealValue("nuTime2D")
                energyTmp = argSet.getRealValue("nuEnergy2D")
                weight = preSelData2D.weight()
                weightError = preSelData2D.weightError()
                newArgSet2D.setRealValue('nuTime2D', timeTmp+dT)
                newArgSet2D.setRealValue('nuEnergy2D', energyTmp)

                fitdata2D.add( newArgSet2D, weight, weightError)
            
            #print('fitdata2D.sumEntries()', round(fitdata2D.sumEntries()))
            res = fitPdf2D.fitTo(fitdata2D, ROOT.RooFit.Range('nllRange'),
                            ROOT.RooFit.Save(), ROOT.RooFit.Extended(),
                            ROOT.RooFit.PrintLevel(-1))

            nllVal = res.minNll()
            #print(nllVal)

            del fitdata2D

            if math.isnan(nllVal) != True:
                tScan.append( dT )
                nllList.append( nllVal )
                nsigList.append(nsig.getVal())
                if locMinNLL>nllVal:
                    locMinNLL = nllVal
                    bestfitDT = dT
        

        print("coarse minimum : ", bestfitDT, locMinNLL )
        print('################ End of Coase Scan ################')

        nllList, tScan, nsigList = [], [], []
        dTexp = bestfitDT    # initial value for deltaT ...
        dTstep = 0.0001
        bestStatus = 0
        NScans = 25
        locMinNLL, bestfitDT, bestfitNSIG = 100.0, 1.0, 0.

        for iScan in range(NScans):
            dT = dTexp + (2*iScan-NScans)*dTstep/2
            if dT>extDTmax or dT<extDTmin:
                print('Warning!: the data may be shifted outside the fitting range')
                continue
            ################## New Data Set
            #newTList, energyList = [], []
            fitdata2D = ROOT.RooDataSet('fitdata2D', 'fitdata2D', newArgSet2D)
            for iEvt in range(nFitEvts2D):
                argSet = preSelData2D.get(iEvt)
                timeTmp = argSet.getRealValue('nuTime2D')
                energyTmp = argSet.getRealValue("nuEnergy")
                weight = preSelData2D.weight()
                weightError = preSelData2D.weightError()
                newArgSet2D.setRealValue('nuTime2D', timeTmp+dT)
                newArgSet2D.setRealValue('nuEnergy2D', energyTmp)

                fitdata2D.add( newArgSet2D, weight, weightError)

                fitdata2D.add( newArgSet2D, weight, weightError)

            nuTime2D.setRange('nllRange', Tmin+dT, Tmax + dT)

            res = fitPdf2D.fitTo(fitdata2D, ROOT.RooFit.Range('nllRange'),
                            ROOT.RooFit.Save(), ROOT.RooFit.Extended(),
                            ROOT.RooFit.PrintLevel(-1))
            nllVal = res.minNll()

            if math.isnan(nllVal) != True:
                tScan.append( dT )
                nllList.append( nllVal )
                nsigList.append(nsig.getVal())
                if locMinNLL>nllVal:
                    locMinNLL = nllVal
                    bestfitDT = dT
                    bestfitNSIG = nsig.getVal()
                    bestStatus = res.status()


        print("Fine minimum : ", bestfitDT, locMinNLL )
        print('################ End of fine Scan ################')


        ## ----- Here is a preliminary checks on fitting quantities ----- #
        
        badFlag = False
        locMinNll, locMaxNll = [], []
        for idd in range(1, len(nllList)-1, 1):
            if nllList[idd-1] < nllList[idd] and nllList[idd] > nllList[idd+1] :
                locMaxNll.append(nllList[idd])
            if nllList[idd-1] > nllList[idd] and nllList[idd] < nllList[idd+1] :
                locMinNll.append(nllList[idd])

        if len(locMaxNll) > 0 or len(locMinNll) > 1 :   # unexpected shape
            badFlag = True
            print("locMinNum = %d, locMaxNum = %d" %(len(locMinNll), len(locMaxNll))) 

        ## -------------------------------------------------------------- #

        del preSelData2D

        for k in range(len(nllList)):
            resfileRaw.write('%d  %8.5f  %12.3f  %12.3f\n'%(iSub, tScan[k], nllList[k], nsigList[k]))

        #resfileSummary.write('%d  %8.5f  %12.3f  %12.3f\n'%(iSub, bestfitDT, locMinNLL, bestfitNSIG))
        resfileSummary.write('%d  %9.6f  %12.3f  %12.3f  %d  %d\n'%(iSub, bestfitDT, locMinNLL, bestfitNSIG, nFitEvts2D, bestStatus))

    resfileRaw.close()
    resfileSummary.close()
    pdfFile.Close()


         #   #if iSub % 5 == 0:
         #    if badFlag == True:
         #       xframe = nuTime.frame(ROOT.RooFit.Title("data group %d, dT: %8.5f"%(iSub, dT)), 
         #                       ROOT.RooFit.Range('nllRange'), ROOT.RooFit.Bins(100))
         #       xframe.GetYaxis().SetTitleOffset(1.4)
         #       can.cd()
         #       fitdata.plotOn(xframe)
         #       ## fitPdf.plotOn(xframe, ROOT.RooFit.ProjectionRange('nllRange'),
         #       ##                 ROOT.RooFit.LineColor(2) )
         #       #                    #ROOT.RooFit.DrawOption("E"),
         #       #                    #ROOT.RooFit.FillColor(ROOT.kOrange),
         #       #                    #ROOT.RooFit.MoveToBack() )
         #       fitPdf1D.plotOn(xframe, ROOT.RooFit.Range('nllRange'),
         #               ROOT.RooFit.Normalization(nFitEvts1D, ROOT.RooAbsReal.NumEvent),
         #               ROOT.RooFit.LineColor(2))
         #       xframe.Draw()
         #       can.Print(pdfilename)

         # fitdata.plotOn(yframe)
         # fitPdf3.plotOn(xframe, ROOT.RooFit.ProjectionRange('fitRange'), 
         #         ROOT.RooFit.Normalization(fitNevts, ROOT.RooAbsReal.NumEvent),
         #         ROOT.RooFit.LineColor(2))


        # Draw the nll profile of iSub toy dataset
        #if badFlag == True:
        ##if iSub % 5 == 0:
        #    tScanX = array('d', tScan)
        #    nllListY = array('d', nllList)
        #    nsigListY = array('d', nsigList)
        #    gr1 = ROOT.TGraph(len(nllList), tScanX, nllListY)
        #    gr2 = ROOT.TGraph(len(nllList), tScanX, nsigListY)
        #    setGraphStyle(gr1, [20, 1.2, 1], [2, 1, 4])
        #    setGraphStyle(gr2, [20, 1.2, 1], [2, 1, 4])
        #    can2.cd(1)
        #    gr1.Draw('APC')
        #    gr1.GetXaxis().SetTitle('#delta_{T} (s)')
        #    gr1.GetYaxis().SetTitle('nll')
        #    can2.cd(2)
        #    gr2.Draw('APC')
        #    gr2.GetXaxis().SetTitle('#delta_{T} (s)')
        #    gr2.GetYaxis().SetTitle('number of events')
        #    can2.Print(pdfilename)




#################################################################################
    #file = open(prefix+'_data%2.1feV_pdf%2.1feV_%2.1fkpc_1DFit.txt'%(nuMass1, nuMass2, dist), 'w')
    #for k in range(len(nllList)):
    #    file.write('%8.5f  %12.3f  %12.3f\n'%(tScan[k], nllList[k], nsigList[k]))
    #file.close()

    #fig2 = plt.figure(figsize=(8,5))
    #ax = fig2.add_subplot()
    #ec = plt.plot( tScan, nllList, 'ro')
    #ax.set_xlabel('$\delta_{T}$ (s)')
    #plt.show()



