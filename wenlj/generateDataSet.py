import ROOT
import math
import matplotlib.pyplot as plt
import sys
import json
import numpy as np
from array import array
from scipy.stats import poisson

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

def printHist(inputHist):
    Tmin, Tmax, Emin, Emax = -0.05, 0.3, 0.2, 5.0
    print('Integral (%4.2f, %4.2f), (%4.2f, %4.2f): %10.3f'
        %(Tmin, Tmax, Emin, Emax, getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")*NEVENTS))

    Tmin, Tmax, Emin, Emax = -0.05, 0.2, 0.2, 3.0
    print('Integral (%4.2f, %4.2f), (%4.2f, %4.2f): %10.3f'
        %(Tmin, Tmax, Emin, Emax, getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")*NEVENTS))
    
    Tmin, Tmax, Emin, Emax = -0.05, 0.3, 0.2, 3.0
    print('Integral (%4.2f, %4.2f), (%4.2f, %4.2f): %10.3f'
        %(Tmin, Tmax, Emin, Emax, getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")*NEVENTS))

if __name__ == "__main__":

    if len(sys.argv) < 8:
        print("Usage: %s [modelNum] [channel] [nuType] [nuMass] [MH] [dist] [Ethr] [group]"  % (sys.argv[0]))
        print("modelNum: check the SNsim/simulation/data directory")
        print("channel : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("nuMass  : neutrino mass that used to generate dataset, in eV unit")
        print("MH      : mass hierarchy, 0 for no osc, 1 for NO, 2 for IO")
        print("dist    : SN distance in kpc unit")
        print("Ethr    : fit Energy threshold")
        print("group    : index of the groups")
        print("example: python generateDataSet.py 82503 1 0 0.0 10.0 0.1 1")
        sys.exit(1)

    ##################################
    # Input variables
    ROOT.RooFit.PrintLevel(-1)
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
    
    nuMass = float( sys.argv[4] )# unit: eV
    print('nuMass for producing dataset : ', nuMass, ' eV')

    MH = int( sys.argv[5] )    
    print("MO :", MH)

    dist = float( sys.argv[6] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[7] )# unit: kpc
    print('fit Ethr: ', Ethr, ' MeV')
    
    group = int( sys.argv[8] )# group number
    print('group number: ', group)

    ##################################
    # Obtain histograms from the root files
    prefix = "./etSpec/fineSpec/TEvisPDF_"
    prefix += 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)
    filename = prefix+'_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25.root'%(nuMass, dist)
    #filename = prefix+'_mNu%2.1feV_%2.1fkpc_0.1s_v4.root'%(nuMass, dist)
    print(filename)
    dataFile  = ROOT.TFile(filename, 'read')
    inputHist = dataFile.Get("hET_mod%d_cha%d_mh%d"%(modelNum, chaname, MH))
    inputBKG = dataFile.Get("hETBKG_mod%d_cha%d_mh%d"%(modelNum, chaname, MH))

    # check the number of integrated events
    print('Integral(), all: ', inputHist.Integral("width"))

    Tmin = inputHist.GetXaxis().GetXmin()
    Tmax = inputHist.GetXaxis().GetXmax()
    Emin = inputHist.GetYaxis().GetXmin()
    Emax = inputHist.GetYaxis().GetXmax()
    print('Tmin, Tmax, Emin, Emax: ', Tmin, Tmax, Emin, Emax)
    nuEnergy = ROOT.RooRealVar("nuEnergy", "nuEnergy", Emin, Emax)
    nuTime   = ROOT.RooRealVar("nuTime", "nuTime", Tmin, Tmax)
    nuEnergy.setRange('fullRange', Emin, Emax)
    nuTime.setRange('fullRange', Tmin, Tmax)
    
    # Integrated Number of Events in the full range of the histograms
    NEVENTS = getIntegral(inputHist, Tmin, Tmax, Emin, Emax ,"width")
    print('NEVENTS in full range: ', NEVENTS)
    NEVENTS = round(NEVENTS)
    print('NEVENTS in full range (round): %10.3f'%(NEVENTS))

    # integrated Number of Events in the ROI
    if MH==0:
        Tmin, Tmax = -0.015, 0.02
    if MH==1:
        #Tmin, Tmax = -0.015, 0.05
        #Tmin, Tmax = -0.04, 0.03
        Tmin, Tmax = -0.04, 0.1
    if MH==2:
        #Tmin, Tmax = -0.015, 0.03
        #Tmin, Tmax = -0.04, 0.03
        Tmin, Tmax = -0.04, 0.1
    subNEVENTS = getIntegral(inputHist, Tmin, Tmax, Ethr, Emax ,"width")
    print('NEVENTS in ROI (Tmin, Tmax, Ethr, Emax): %f in (%5.2f, %5.2f, %5.2f, %5.2f)'%\
          (subNEVENTS, Tmin, Tmax, Ethr, Emax))
    Emax = 10.
    subNEVENTS = getIntegral(inputHist, Tmin, Tmax, Ethr, Emax ,"width")
    print('NEVENTS in ROI (Tmin, Tmax, Ethr, Emax): %f in (%5.2f, %5.2f, %5.2f, %5.2f)'%\
          (subNEVENTS, Tmin, Tmax, Ethr, Emax))
    subNEVENTS = round(subNEVENTS)
    print('NEVENTS in ROI (round): ', subNEVENTS)

    # 1D histogram in time domain
    binEthr = inputHist.GetYaxis().FindBin(Ethr)
    print('Ethr, binEthr, lower bin edge of binEthr: ', 
          Ethr, ',', binEthr, ', ', inputHist.GetYaxis().GetBinLowEdge(binEthr)) 

    binEvismax = inputHist.GetYaxis().FindBin(Emax)
    print('Ethr, binEvismax, lower bin edge of binEvismax: ', 
          Ethr, ',', binEvismax, ', ', inputHist.GetYaxis().GetBinLowEdge(binEvismax)) 

    inputHist1D = inputHist.ProjectionX('_pdfsigpx', binEthr, binEvismax, '')
    inputBKG1D  =  inputBKG.ProjectionX('_pdfbkgpx', binEthr, binEvismax, '')

    #NEVENTS = NEVENTS * 4 # special test for 10 kpc, group # == 999
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))
    #subNEVENTS = subNEVENTS * 4 # special test for 10 kpc, group # == 999
    #print('Modified NEVENTS in ROI (round): ', subNEVENTS)

    #NEVENTS = NEVENTS * 25 # special test for 10 kpc, group # == 88X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))
    #subNEVENTS = subNEVENTS * 25 # special test for 10 kpc, group # == 88X
    #print('Modified NEVENTS in ROI (round): ', subNEVENTS)

    #NEVENTS = NEVENTS * 100 # special test for 10 kpc, group # == 77X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))
    #subNEVENTS = subNEVENTS * 100 # special test for 10 kpc, group # == 77X
    #print('Modified NEVENTS in ROI (round): ', subNEVENTS)

    #NEVENTS = int(NEVENTS * 6.25) # special test for 5 kpc, group # == 88X
    #print('Modified NEVENTS in full range (round): %10.3f'%(NEVENTS))
    #subNEVENTS = int(subNEVENTS * 6.25) # special test for 5 kpc, group # == 88X
    #print('Modified NEVENTS in ROI (round): ', subNEVENTS)

    # Poisson Statistics
    # nevtList = []
    # # rvs(mu, loc=0, size=1, random_state=None)
    # for i in range(100):
    #     sampleNevt = poisson.rvs(subNEVENTS)
    #     nevtList.append( sampleNevt )
    nevtList = poisson.rvs(subNEVENTS, size=500, random_state=group)
    nevtStats = len(nevtList)
    print(nevtList)

    # normFactor = 1./NEVENTS
    # inputHist.Scale(normFactor)
    # print('normalized all: ', getIntegral(inputHist, Tmin, Tmax, Emin, Emax, "width"))

    ################# Generate PDFs
    # generate RooDataHist
    dataHist = ROOT.RooDataHist('dataHist', 'dataHist', 
                ROOT.RooArgList(nuTime, nuEnergy), 
                ROOT.RooFit.Import(inputHist, True))
    
    dataHist1D = ROOT.RooDataHist('dataHist1D', 'dataHist1D', 
                ROOT.RooArgList(nuTime), 
                ROOT.RooFit.Import(inputHist1D, True))

    # Obtain the data PDF
    dataPdf = ROOT.RooHistPdf("dataPdf", "dataPdf", 
                ROOT.RooArgSet(nuTime, nuEnergy), dataHist, 2)
    
    dataPdf1D = ROOT.RooHistPdf("dataPdf1D", "dataPdf1D", 
                ROOT.RooArgSet(nuTime), dataHist1D, 2)

    ################# Generate data set from the data PDF
    ROOT.RooRandom.randomGenerator().SetSeed(group)

    #data = dataPdf.generate(ROOT.RooArgSet(nuTime, nuEnergy), nevtStats*NEVENTS)
    #data1D = dataPdf1D.generate(ROOT.RooArgSet(nuTime), nevtStats*NEVENTS)
    #
    #print('data.sumEntries(): ', round(data.sumEntries()) )
    #print('data1D.sumEntries(): ', round(data1D.sumEntries()) )

    # Modify the variable range
    nuEnergy.setRange('preSelRange', Ethr, Emax)
    nuTime.setRange('preSelRange', Tmin, Tmax)

    #preSelData1D = data1D.reduce(ROOT.RooFit.CutRange('preSelRange'))
    #nFitEvts1D = round(preSelData1D.sumEntries())
    #print('preSelData1D.sumEntries(): ', nFitEvts1D )

    #preSelData = data.reduce(ROOT.RooFit.CutRange('preSelRange'))
    #nFitEvts = round(preSelData.sumEntries())
    #print('preSelData.sumEntries(): ', nFitEvts )

    # ################# Save dataset to ROOT file
    outPrefix = "./dataset/fineSpec/%2.1fkpc/TEvisDATA_"%(dist)
    outPrefix += 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)
    filename = outPrefix+'_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25_Ethr%2.1fMeV_group%d.root'%(nuMass, dist, Ethr, group)
    outFile = ROOT.TFile(filename, 'recreate')
    # tinfo = ROOT.TTree("evtTree", "data set")
    # tinfo.Branch("model", imod, "model/I")
    # tinfo.Branch("channel", icha, "channel/I")
    # tinfo.Branch("tStart", tStart, "tStart/D")
    # tinfo.Branch("tEnd", tEnd, "tEnd/D")
    
    # ################# Fill the tree with fixed statistics
    timeF1D = array('d', [0.])
    timeF   = array('d', [0.])
    energyF = array('d', [0.])
    evtID = array('i', [0])
    tFixedStat = ROOT.TTree("tFixedStat", "data set")
    tFixedStat.Branch("nuTime1D", timeF1D, "nuTime1D/D")
    tFixedStat.Branch("nuTime2D", timeF, "nuTime2D/D")
    tFixedStat.Branch("nuEnergy", energyF, "nuEnergy/D")
    tFixedStat.Branch("evtID", evtID, "evtID/I")

    # fixed statistics
    iEvt = 0
    tListTest = []
    for i in range(nevtStats):
        if i%100==0:
            print('nevtStats', i)

        data = dataPdf.generate(ROOT.RooArgSet(nuTime, nuEnergy), 2*NEVENTS)
        data1D = dataPdf1D.generate(ROOT.RooArgSet(nuTime), 2*NEVENTS)
        preSelData1D = data1D.reduce(ROOT.RooFit.CutRange('preSelRange'))
        preSelData = data.reduce(ROOT.RooFit.CutRange('preSelRange'))

        evtID[0] = i + 1
        for k in range(subNEVENTS):
            argSet  = preSelData1D.get(k)
            timeF1D[0]   = argSet.getRealValue('nuTime')
            #if k<20:
            #    tListTest.append( timeF1D[0] )

            argSet  = preSelData.get(k)
            timeF[0]   = argSet.getRealValue('nuTime')
            energyF[0] = argSet.getRealValue('nuEnergy')
 
            tFixedStat.Fill()

        del data
        del data1D
        del preSelData1D
        del preSelData
    tFixedStat.Write()
    #print(tListTest)

    ################# Fill the tree with Poisson statistics
    timeR1D = array('d', [0.])
    timeR   = array('d', [0.])
    energyR = array('d', [0.])
    tRndStat = ROOT.TTree("tRndStat", "data set")
    tRndStat.Branch("nuTime1D", timeR1D, "nuTime1D/D")
    tRndStat.Branch("nuTime2D", timeR, "nuTime2D/D")
    tRndStat.Branch("nuEnergy", energyR, "nuEnergy/D")
    tRndStat.Branch("evtID", evtID, "evtID/I")

    iEvt = 0
    for i in range(nevtStats):
        if i%100==0:
            print('nevtStats', i)

        data = dataPdf.generate(ROOT.RooArgSet(nuTime, nuEnergy), 2*NEVENTS)
        data1D = dataPdf1D.generate(ROOT.RooArgSet(nuTime), 2*NEVENTS)
        preSelData1D = data1D.reduce(ROOT.RooFit.CutRange('preSelRange'))
        preSelData = data.reduce(ROOT.RooFit.CutRange('preSelRange'))

        evtID[0] = i + 1
        nEvt = nevtList[i]
        for k in range(nEvt):
            argSet  = preSelData1D.get(k)
            timeR1D[0]   = argSet.getRealValue('nuTime')

            argSet  = preSelData.get(k)
            timeR[0]   = argSet.getRealValue('nuTime')
            energyR[0] = argSet.getRealValue('nuEnergy')
            tRndStat.Fill()

        del data
        del data1D
        del preSelData1D
        del preSelData
    tRndStat.Write()

    outFile.Close()
