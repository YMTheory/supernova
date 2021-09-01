### from Wikipedia
### The two-sample K–S test is one of the most useful and general 
### nonparametric methods for comparing two samples, as it is 
### sensitive to differences in both location and shape of the 
### empirical cumulative distribution functions of the two samples

import numpy as np
import matplotlib.pyplot as plt
import ROOT
from scipy import stats
import sys, os
import namespace as ns
import uproot as up
import boost_histogram as bh


# empirical distribution function Fn
def edf(x, samples, deltaT):
    Fn = 0.
    for i in samples :
        if i + deltaT <= x:
            Fn += 1
    return Fn/len(samples)


# cumulative distribution function Fx :
def cdf(x, hist, xmin, xmax):
    bin1 = hist.FindBin(xmin)
    bin2 = hist.FindBin(xmax)

    bb = hist.FindBin(x)

    integ_total = hist.Integral(bin1, bin2)
    integ_part  = hist.Integral(bin1, bb)

    return integ_part / integ_total


def getValue(strline):
    values = strline.split()
    x1 = int(values[0])
    x2 = float(values[1])
    x3 = float(values[2])
    x4 = float(values[3])
    return (x1, x2, x3, x4)


def getNLL(filename, num_datasets):
    groupNum, deltaT, nllVal, nsig = [], [], [], []
    
    if os.path.exists(filename) == False:
        print("%s does not exist !" %filename)
        return groupNum, deltaT, nllVal, nsig

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()

    if num_datasets > len(lines):
        print(filename)
        print('Error! num_datasets is larger than the number of rows in txt file!')
        return groupNum, deltaT, nllVal, nsig

    for k in range(num_datasets):
        (num, dT, val, nEvt) = getValue( lines[k] )
        groupNum.append(num)
        deltaT.append(dT)
        nllVal.append(val)
        nsig.append(nEvt)
    return groupNum, deltaT, nllVal, nsig






if __name__ == "__main__" :
    if len(sys.argv) < 11:
        print("Usage: %s [modelNum] [channel] [nuType] [dataNuMass] [dataMH] [fitNuMass] [pdfMH] [dist] [Ethr] [group]"%(sys.argv[0]))
        print("modelNum: check the SNsim/simulation/data directory")
        print("channel : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("dataNuMass : neutrino mass that used to generate dataset, in eV unit")
        print("dataMH  : NMO that used to generate dataset, 0 for NoOsc, 1 for NH, 2 for IH")
        print("fitNuMass : neutrino mass that used in PDF to fit the dataset, in eV unit")
        print("pdfMH   : NMO that used to generate PDF")
        print("dist    : SN distance in kpc unit")
        print("Ethr    : fit Energy threshold")
        print("group   : index of the groups")
        print("example : python nllFit.py 82503 1 0 0.0 2 0.2 1 10.0 0.1 1")
        sys.exit(1)

    ##################################
    # choose SN model
    #modelNum = 82503
    modelNum = int( sys.argv[1] )

    chaName = ['nup','nue','IBD']
    chaname = int( sys.argv[2] )
    print("channelName: ", chaname, chaName[chaname])

    nuType = int( sys.argv[3] )
    # -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x
    nuTypeName = ['nu_e', 'anti_nue', 'nu_x', 'all']
    print("neutrino type: ", nuType, nuTypeName[nuType])
    
    dataNuMass = float( sys.argv[4] )# unit: eV
    print('nuMass for producing dataset : ', dataNuMass, ' eV')

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


    fitTmin, fitTmax = -0.015, 0.02

    # dataset file loading 
    datafilename = ns.dataFileName(modelNum, dist, chaname, nuType, dataNuMass, dataMH, Ethr, group)
    inputTree = up.open(datafilename)["tFixedStat"]
    evtId   = inputTree["evtID"].array()
    dataset = inputTree["nuTime1D"].array()

    for iEvt in range(1, 11, 1):
        print(" -----> Processing event %d" %iEvt)

        nuTimeArr = []
        for i, j in zip(evtId, dataset) :
            if i == iEvt and j >= fitTmin and j<= fitTmax:
                nuTimeArr.append(j)

        nuMassList = np.arange(0., 2.0, 0.1)
        Dmax = np.zeros(20)

        iPDF = 0
        for fitNuMass in nuMassList:
            # pdf file loading 
            pdffilename = ns.pdfSumFileName(modelNum, dist, chaname, fitNuMass, fitMH, nuType, 25.00)
            pdfFile  = ROOT.TFile(pdffilename, 'read')
            inputHist = pdfFile.Get("TEvisPdf")
            # 1D histogram in time domain
            binEthr = inputHist.GetYaxis().FindBin(Ethr)
            binEmax = inputHist.GetYaxis().FindBin(4.)
            ################# Input 1D Signal/Background 1D Histograms
            inputHist1D = inputHist.ProjectionX('_pdfsigpx', binEthr, binEmax, '')
            

            # fitting results :
            resfileSummary = ns.fitResFileName(modelNum, chaname, dataMH, dataNuMass, fitMH, fitNuMass, dist, Ethr, group)[1]

            if not os.path.exists(resfileSummary) :
                print("%s does not exists ! " %resfileSummary)
                sys.exit(1)
            

            nStatsPerGroup = 100

            groupNum, deltaT, nllVal, nsig = [], [], [], []
            t1, t2, t3, t4 = getNLL(resfileSummary, nStatsPerGroup)
            groupNum = groupNum + t1
            deltaT = deltaT + t2
            nllVal = nllVal + t3
            nsig   = nsig   + t4

            dT = deltaT[iEvt]

            xList = np.arange(fitTmin+dT, fitTmax+dT, 0.00005)
            yList = []
            for i in xList:
                yList.append( cdf(i, inputHist1D, fitTmin+dT, fitTmax+dT) )

            dList = []
            for i in xList :
                dList.append(edf(i, nuTimeArr, dT))

            yList = np.array(yList)
            dList = np.array(dList)
            deltaCDF = np.abs(yList - dList)
            
            Dmax[iPDF] = np.max(deltaCDF)
            iPDF += 1


        plt.plot(nuMassList, Dmax, "o-")
    plt.show()

    #plt.plot(xList, yList, color="blue", label="model CDF")
    #plt.plot(xList, dList, color="red",  label="empirical CDF")

    #plt.xlabel("time/s")
    #plt.ylabel("CDF")
    #plt.show()

    









