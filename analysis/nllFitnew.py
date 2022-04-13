import ROOT
import math
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import sys
import json
import scipy
import numpy as np
from array import array

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

    if len(sys.argv) < 12:
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
        print("fitEmax : fit Maximum Energy")
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

    fitEmax = float( sys.argv[10] )# unit: kpc
    print('fit fitEmax: ', fitEmax, ' MeV')

    group = int( sys.argv[11] )# group number
    print('group number: ', group)


    # ###################################
    # # Obtain the fitting PDF
    #filename = ns.pdfSumFileName(modelNum, 10.0, chaname, nuMass2, fitMH, nuType, 25.00)
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/5kpc/evistSpec_mod%d_cha%d_nuMass%.1f.root" %(modelNum, chaname, nuMass2)
    #filename = ns.pdfSumFileName(modelNum, dist, chaname, nuMass2, fitMH, nuType, 25.00)
    print(filename)
    pdfFile  = ROOT.TFile(filename, 'read')
    #inputHist = pdfFile.Get("TEvisPdf")
    if chaname == 1:
        inputHist = pdfFile.Get("hET_mod82503_cha%d_mh%d"%(chaname, fitMH))
    #inputBKG  = pdfFile.Get("TBkgEvisPdf")
    if chaname == 0:
        inputHist = pdfFile.Get("hET_mod82503_cha%d_mh%d"%(chaname, 1))


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

    #NEVENTS = NEVENTS * 25 # special test for 10.0pc, group # == 88X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    #NEVENTS = NEVENTS * 100 # special test for 10 kpc, group # == 77X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))

    # 1D histogram in time domain
    binEthr = inputHist.GetYaxis().FindBin(Ethr)
    print('Ethr, binEthr, lower bin edge of binEthr: ', 
          Ethr, ',', binEthr, ', ', inputHist.GetYaxis().GetBinLowEdge(binEthr)) 
    binEmax = inputHist.GetYaxis().FindBin(fitEmax)
    print('fitEmax, binEmax, lower bin edge of binEmax: ', 
          fitEmax, ',', binEmax, ', ', inputHist.GetYaxis().GetBinLowEdge(binEmax)) 

    ################# Input 1D Signal/Background 1D Histograms

    inputHist1D = inputHist.ProjectionX('_pdfsigpx', binEthr, binEmax, '')
    #inputBKG1D  =  inputBKG.ProjectionX('_pdfbkgpx', binEthr, -1, '')

    sigHist = ROOT.RooDataHist('sigHist', 'sigHist', 
          ROOT.RooArgList(nuTime, nuEnergy), 
          ROOT.RooFit.Import(inputHist, True) )

    #bkgHist = ROOT.RooDataHist('bkgHist', 'bkgHist', 
    #      ROOT.RooArgList(nuTime, nuEnergy), 
    #      ROOT.RooFit.Import(inputBKG, True) )

    sigPdf = ROOT.RooHistPdf('sigPdf', 'sigPdf', 
          ROOT.RooArgList(nuTime, nuEnergy), sigHist, 2)

    #bkgPdf = ROOT.RooHistPdf('bkgPdf', 'bkgPdf', 
    #      ROOT.RooArgList(nuTime, nuEnergy), bkgHist, 2)

    # 1D PDF in time domain
    sigHist1D = ROOT.RooDataHist('sigHist1D', 'sigHist1D', 
          ROOT.RooArgList(nuTime), 
          ROOT.RooFit.Import(inputHist1D, True) )

    #bkgHist1D = ROOT.RooDataHist('bkgHist1D', 'bkgHist1D', 
    #      ROOT.RooArgList(nuTime), 
    #      ROOT.RooFit.Import(inputBKG1D, True) )

    sigPdf1D = ROOT.RooHistPdf('sigPdf1D', 'sigPdf1D', 
          ROOT.RooArgList(nuTime), sigHist1D, 2)

    #bkgPdf1D = ROOT.RooHistPdf('bkgPdf1D', 'bkgPdf1D', 
    #      ROOT.RooArgList(nuTime), bkgHist1D, 2)

    #scaleFac = ROOT.RooRealVar('scaleFac', 'scaleFac', )
    nsig = ROOT.RooRealVar("nsig","#signal events", NEVENTS, 0., NEVENTS*10000)
    #nbkg = ROOT.RooRealVar("nbkg","#background events", 0.0)

    #fitPdf = ROOT.RooAddPdf('fitPdf', 'fitPdf', 
    #      ROOT.RooArgList(sigPdf, bkgPdf), 
    #      ROOT.RooArgList(nsig, nbkg) )

    # fitPdf1D from nllFit.py
    fitPdf1D = ROOT.RooExtendPdf("fitPdf1D", "fitPdf1D", sigPdf1D, nsig)
    print("1DfitPdf has been constructed successfully !")

    # old fitPdf1D in this script:
    #fitPdf1D = ROOT.RooAddPdf('fitPdf1D', 'fitPdf1D', 
    #      ROOT.RooArgList(sigPdf1D, bkgPdf1D), 
    #      ROOT.RooArgList(nsig, nbkg) )

    ##################################
    # Obtain data sets from the root files
    if dataMH == 1:
        dmo = "NO"
    if dataMH == 2:
        dmo = "IO"
    #filename = "/junofs/users/miaoyu/supernova/miao/mini-pkg/fitter/nuMass0.0eV_"+dmo+"_stat10000_Emax10MeV.root"
    #filename = ns.dataFileName(modelNum, dist, chaname, nuType, nuMass1, dataMH, Ethr, group)
    #filename = ns.dataOldFileName(modelNum, dist, chaname, nuMass1, dataMH, Ethr, group)
    filename = "/junofs/users/miaoyu/supernova/production/Data/5kpc/data_mod%d_cha%d_mh%d_nuMass%.1f.root"%(modelNum, chaname, dataMH, nuMass1)
    print(filename)
    dataFile  = ROOT.TFile(filename, 'read')
    inputTree = dataFile.Get('evt')

    # change the range for pre-selection to the dataset
    if chaname == 2:  # snDet.IBD
        Tmin, Tmax, Emin, Emax = -0.02, 0.2, 1.0, 10.0
    if chaname == 1:  # snDet.NuE
        #Tmin, Tmax, Emin, Emax = -0.02, 0.025, 0.1, 4.0
        Tmin, Tmax, Emin, Emax = -0.015, 0.02, Ethr, fitEmax
    if chaname == 0:
        Tmin, Tmax, Emin, Emax = -0.015, 0.02, Ethr, fitEmax

    nuTime1D = ROOT.RooRealVar("nuTime1D", "nuTime1D", Tmin, Tmax)
    nuEnergy = ROOT.RooRealVar("nuEnergy", "nuEnergy", Emin, Emax)
    evtID    = ROOT.RooRealVar("evtID", "evtID", 0, 500)

    resPath = "/junofs/users/miaoyu/supernova/analysis/fitResults/"
    resFiles = ["res_mod%d_cha%d_mh%d%d_nuMass%.1feV_stat10000.txt"%(modelNum, chaname, dataMH, fitMH, nuMass2), "res_mod%d_cha%d_mh%d%d_nuMass%.1feV_stat10000_summary.txt"%(modelNum, chaname, dataMH, fitMH, nuMass2)]
    #resFiles = ["test.txt", "summary.txt"]
    resfileRaw = open(resPath+resFiles[0], "w")
    resfileSummary = open(resPath+resFiles[1], 'w')

    can = ROOT.TCanvas("can", "can", 1200, 800)
    pdfilename = resPath + "res_mod%d_cha%d_mh%d%d_nuMass%.1feV_stat10000.pdf"%(modelNum, chaname, dataMH, fitMH, nuMass2)
    #pdfilename = resPath + "draw.pdf"
    can.Print(pdfilename + '[')

    can2 = ROOT.TCanvas("can2", "can2", 1200, 800)
    can2.Divide(2,1)
    can2.cd(1).SetLeftMargin(0.15)
    can2.cd(1).SetRightMargin(0.05)
    can2.cd(2).SetRightMargin(0.05)

    newArgSet = ROOT.RooArgSet(nuTime)


    nStat = 500
    for iSub in range(0, nStat):
        Tmin, Tmax = -0.015, 0.02     # Fitting time range
    #for iSub in range(29, 30):
        print('Process toy dataset %d'%iSub)
        # Get RooDataSet from the tree
        preSelData1D = ROOT.RooDataSet("dataGroup%d"%(iSub), "dataset with (e,t)", 
                        inputTree, ROOT.RooArgSet(nuTime1D, evtID), 
                        "evtID==%d"%(iSub) )

    
        nFitEvts1D = round(preSelData1D.sumEntries())
        print('preSelData1D.sumEntries(): ', nFitEvts1D )

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

        for iScan in range(NScans):
            dT = dTexp + (2*iScan-NScans)*dTstep/2
            if dT>extDTmax or dT<extDTmin:
                print('Warning!: the data may be shifted outside the fitting range')
                continue

            ################## New Data Set
            #print('################ iSub %d, Start of Coase Scan %d ################'%(iSub, iScan))
            fitdata = ROOT.RooDataSet('fitdata', 'fitdata', newArgSet)
            #print('fitdata.sumEntries()', round(fitdata.sumEntries()))
            # change the range for NLL fitting the dataset
            nuTime.setRange('nllRange', Tmin + dT, Tmax + dT)
            nuEnergy.setRange('nllRange', Emin, Emax)
        
            #newTList, energyList = [], []
            #nNu_nllRange = 0
            for iEvt in range(nFitEvts1D):
                argSet = preSelData1D.get(iEvt)
                timeTmp = argSet.getRealValue('nuTime1D')
                #if timeTmp >-0.015 and timeTmp <= 0.02:
                #    nNu_nllRange += 1
                #energyTmp = argSet.getRealValue('nuEnergy')
                weight = preSelData1D.weight()
                weightError = preSelData1D.weightError()
                #print('before, time, energy, weight, weightError: ', 
                #        timeTmp, energyTmp, weight, weightError)

                newArgSet.setRealValue('nuTime', timeTmp + dT)
                #newArgSet.setRealValue('nuEnergy', energyTmp)
                #newTList.append(timeTmp + dT)
                #energyList.append(energyTmp)
                        
                #print('change, timeTmp: ', newArgSet.getRealValue('nuTime1D'))
                fitdata.add( newArgSet, weight, weightError)

            #print("Neutrino number in ROI : %d"%nNu_nllRange)
            res = fitPdf1D.fitTo(fitdata, ROOT.RooFit.Range('nllRange'),
                            ROOT.RooFit.Save(), ROOT.RooFit.Extended(),
                            ROOT.RooFit.PrintLevel(-1))
            nllVal = res.minNll()
            
            #print("final sat. neutrino number : %d" %nNu_nllRange)
            #print('fitdata.sumEntries()', round(fitdata.sumEntries()))
            # print('################ iSub %d, End of Coase Scan %d ################'%(iSub, iScan))
            del fitdata

            if math.isnan(nllVal) != True:
                tScan.append( dT )
                nllList.append( nllVal )
                nsigList.append(nsig.getVal())
                if locMinNLL>nllVal:
                    locMinNLL = nllVal
                    bestfitDT = dT
        # print('coarse minimum', bestfitDT, locMinNLL)
        # print('################ End of Coase Scan ################')
        ###################################################################################
        ################# Likelihood Fine scan
        nllList, tScan, nsigList = [], [], []
        dTexp = bestfitDT
        print('dTexp: ', dTexp)
        dTstep = 0.0001
        NScans = 25
        locMinNLL, bestfitDT, bestfitNSIG = 100.0, 1.0, 0.
        bestStatus = 0

        for iScan in range(NScans):   
            dT = dTexp + (2*iScan-NScans)*dTstep/2
            if dT>extDTmax or dT<extDTmin:
                print('Warning!: the data may be shifted outside the fitting range')
                continue
        
            ################## New Data Set
            #newTList, energyList = [], []
            fitdata = ROOT.RooDataSet('fitdata', 'fitdata', newArgSet)
            for iEvt in range(nFitEvts1D):
                argSet = preSelData1D.get(iEvt)
                timeTmp = argSet.getRealValue('nuTime1D')
                #energyTmp = argSet.getRealValue('nuEnergy')
                weight = preSelData1D.weight()
                weightError = preSelData1D.weightError()

                newArgSet.setRealValue('nuTime', timeTmp + dT)
                #newArgSet.setRealValue('nuEnergy', energyTmp)
                #newTList.append(timeTmp + dT)
                #energyList.append(energyTmp)
                        
                fitdata.add( newArgSet, weight, weightError)

            nuTime.setRange('nllRange', Tmin+dT, Tmax + dT)

            res = fitPdf1D.fitTo(fitdata, ROOT.RooFit.Range('nllRange'),
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


            xframe = nuTime.frame(ROOT.RooFit.Title("data group %d, dT: %8.5f"%(iSub, dT)), 
                            ROOT.RooFit.Range('nllRange'), ROOT.RooFit.Bins(100))
            xframe.GetYaxis().SetTitleOffset(1.4)
            can.cd()
            fitdata.plotOn(xframe)
            fitPdf1D.plotOn(xframe, ROOT.RooFit.Range('nllRange'),
                     ROOT.RooFit.Normalization(nFitEvts1D, ROOT.RooAbsReal.NumEvent),
                     ROOT.RooFit.LineColor(2))
            xframe.Draw()
            can.Print(pdfilename)
            del fitdata

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




         #   #if iSub % 5 == 0:
         #    if badFlag == True:
        #xframe = nuTime.frame(ROOT.RooFit.Title("data group %d, dT: %8.5f"%(iSub, dT)), 
        #                ROOT.RooFit.Range('nllRange'), ROOT.RooFit.Bins(100))
        #xframe.GetYaxis().SetTitleOffset(1.4)
        #can.cd()
        #fitdata.plotOn(xframe)
        #fitPdf1D.plotOn(xframe, ROOT.RooFit.Range('nllRange'),
        #         ROOT.RooFit.Normalization(nFitEvts1D, ROOT.RooAbsReal.NumEvent),
        #         ROOT.RooFit.LineColor(2))
        #xframe.Draw()
        #can.Print(pdfilename)

        #del fitdata


        #print('nll: ', nllList)
        #print(preSelData1D)
        print('################ End of Fine Scan ################')
        del preSelData1D    


        # Draw the nll profile of iSub toy dataset
        #if badFlag == True:
        if iSub % 1 == 0:
            tScanX = array('d', tScan)
            nllListY = array('d', nllList)
            nsigListY = array('d', nsigList)
            gr1 = ROOT.TGraph(len(nllList), tScanX, nllListY)
            gr2 = ROOT.TGraph(len(nllList), tScanX, nsigListY)
            setGraphStyle(gr1, [20, 1.2, 1], [2, 1, 4])
            setGraphStyle(gr2, [20, 1.2, 1], [2, 1, 4])
            can2.cd(1)
            gr1.Draw('APC')
            gr1.GetXaxis().SetTitle('#delta_{T} (s)')
            gr1.GetYaxis().SetTitle('nll')
            can2.cd(2)
            gr2.Draw('APC')
            gr2.GetXaxis().SetTitle('#delta_{T} (s)')
            gr2.GetYaxis().SetTitle('number of events')
            can2.Print(pdfilename)

        # save the results
        for k in range(len(nllList)):
            resfileRaw.write('%d  %8.5f  %12.3f  %12.3f\n'%(iSub, tScan[k], nllList[k], nsigList[k]))

        #resfileSummary.write('%d  %8.5f  %12.3f  %12.3f\n'%(iSub, bestfitDT, locMinNLL, bestfitNSIG))
        resfileSummary.write('%d  %9.6f  %12.3f  %12.3f  %d  %d\n'%(iSub, bestfitDT, locMinNLL, bestfitNSIG, nFitEvts1D, bestStatus))

    resfileRaw.close()
    resfileSummary.close()
    pdfFile.Close()
    can.Print(pdfilename+']')


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
