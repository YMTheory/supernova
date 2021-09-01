import numpy as np
import matplotlib.pyplot as plt
import namespace as ns
import ROOT
import sys, os
import uproot as up

def readBestFit(filename, iEvt, nuMass):
    ff = ROOT.TFile(filename, "read")
    g1name = "Evt" + str(iEvt) + "nuMass" + str(int(nuMass * 10)) + "NLL"
    g1 = ff.Get(g1name)
    g2name = "Evt" + str(iEvt) + "nuMass" + str(int(nuMass * 10)) + "Nsig"
    g2 = ff.Get(g2name)
    
    deltaT, chi2, nsig = [], [], []
    for i in range(g1.GetN()):
        deltaT.append(g1.GetPointX(i))
        chi2.append(g1.GetPointY(i))
        nsig.append(g2.GetPointY(i))

    deltaT = np.array(deltaT)
    chi2 = np.array(chi2)
    nsig = np.array(nsig)

    return deltaT[np.argmin(chi2)], nsig[np.argmin(chi2)], np.min(chi2)



if __name__ == "__main__" :

    if len(sys.argv) < 12:
        print("Usage: %s [modelNum] [channel] [nuType] [nuMass1] [dataMH] [nuMass2low] [nuMass2high] [pdfMH] [dist] [Ethr] [group]"%(sys.argv[0]))
        print("modelNum    : check the SNsim/simulation/data directory")
        print("channel     : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType      : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("nuMass1     : neutrino mass that used to generate dataset, in eV unit")
        print("dataMH      : NMO that used to generate dataset, 0 for NoOsc, 1 for NH, 2 for IH")
        print("nuMass2low  : low neutrino mass that used in PDF to fit the dataset, in eV unit")
        print("nuMass2high : high neutrino mass that used in PDF to fit the dataset, in eV unit")
        print("pdfMH       : NMO that used to generate PDF")
        print("dist        : SN distance in kpc unit")
        print("Ethr        : fit Energy threshold")
        print("group       : index of the groups")
        print("example     : python nllFit.py 82503 1 0 0.0 2 0.2 1 10.0 0.1 1")
        sys.exit(1)

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

    nuMass2low = float( sys.argv[6] )# unit: eV
    print('low nuMass for producing PDF: ', nuMass2low, ' eV')

    nuMass2high = float( sys.argv[7] )# unit: eV
    print('high nuMass for producing PDF: ', nuMass2high, ' eV')

    fitMH = int( sys.argv[8] )
    print('NMO for producing PDF: ', fitMH, MHName[fitMH])

    dist = float( sys.argv[9] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[10] )# unit: kpc
    print('fit Ethr: ', Ethr, ' MeV')

    group = int( sys.argv[11] )# group number
    print('group number: ', group)

    if group == 1:
        statScale = 1.
    if group == 999:
        statScale = 4
    if group == 880:
        statScale = 25
    if group == 770:
        statScale = 100
   
    iEvt = 1

    ##################################
    # Obtain data sets from the root files
    #filename = ns.dataFileName(modelNum, dist, chaname, nuType, nuMass1, dataMH, Ethr, group)
    #filename = ns.dataOldFileName(modelNum, dist, chaname, nuMass1, dataMH, Ethr, group)
    filename = "/junofs/users/miaoyu/supernova/miao/mini-pkg/fitter/nuMass0.0eV_NO_stat1000.root"
    print(filename)
    inputTree = up.open(filename)["tFixedStat"]
    nuTime1D = inputTree["nuTime1D"].array()
    evtId    = inputTree["evtID"].array()

    nuTimeArr = []
    for i, j in zip(evtId, nuTime1D) :
        if i == iEvt :
            nuTimeArr.append(j)



    ##################################
    # Obtain PDF from root files :

    Emin, Emax = Ethr, 4.0
    nPDF = int((nuMass2high - nuMass2low))
    print("Total fitting Pdf number = %d" %nPDF)
    TPdf = [ROOT.TH1D("Tdf%d"%i, "", 4000, -0.1, 0.1) for i in range(nPDF)]
    for massid in range(nPDF):
        nuMass2 = massid / 10. + nuMass2low / 10.
        print("PDF mass now : %.1feV -> Id%d" %(nuMass2, massid))
        filename = ns.pdfSumFileName(modelNum, dist, chaname, nuMass2, fitMH, nuType, 25.00)
        #filename = "./etSpec/fineSpec/TEvisPDF_mod82503_cha1nue_mh2_mNu%.1feV_10.0kpc_0.1s_Evmax25.root"%(nuMass2)
        print(filename)
        pdfFile  = ROOT.TFile(filename, 'read')
        inputHist = pdfFile.Get("TEvisPdf")
        #inputHist = pdfFile.Get("hET_mod82503_cha1_mh2")
        
        Tmin, Tmax = -0.015 , 0.020
        x1, x2 = inputHist.GetXaxis().FindBin(Tmin), inputHist.GetXaxis().FindBin(Tmax)
        y1, y2 = inputHist.GetYaxis().FindBin(Emin), inputHist.GetYaxis().FindBin(Emax)
        integROI = inputHist.Integral(x1, x2, y1, y2, "width")
        statScale = 1000 / integROI
        print(">>>>>>>>>>>>>> statScale : " , statScale)
        
        # 1D projection manually
        for i in range(inputHist.GetNbinsX()):
            contSum = 0.
            for j in range(inputHist.GetNbinsY()):
                tmpE = inputHist.GetYaxis().GetBinCenter(j+1)
                if tmpE < Ethr or tmpE > Emax :
                    continue
                cont = inputHist.GetBinContent(i+1, j+1) * inputHist.GetYaxis().GetBinWidth(j+1)
                contSum += cont
            contSum *= statScale
            TPdf[massid].SetBinContent(i+1, contSum)


    #### return the best fit value :
    #for nuMass2 in np.arange(nuMass2low, nuMass2high, 0.1):
    #    bestFitT, bestFitChi2, bestFitNsig = readBestFit("/junofs/users/miaoyu/supernova/miao/mini-pkg/fitter/check.root", iEvt, nuMass2)
    #    print(nuMass2, bestFitT, bestFitChi2, bestFitNsig)


    #### Draw Fitting Results from TMinuit...
    iMass = 0
    bestFitT = [0, 6.34568e-5, 0.000153827, 0.000371111]
    for iPDF in [0]:
    #for iPDF in [0, 1, 2, 3]:
        nBins_data = 120
        hData = ROOT.TH1D("dataset", "", nBins_data, -0.04, 0.04)
        bin_width = 0.08 / nBins_data
        fitTmin, fitTmax = -0.04, 0.04
        for nuT in nuTimeArr:
            if nuT < fitTmin or nuT > fitTmax :
                continue
            nuT += bestFitT[iMass]
            hData.Fill(nuT)

        print("Total signal number in ROI : %d" %hData.GetEntries())
        data_center, data_content = [], []
        for i in range(nBins_data):
            if hData.GetBinContent(i+1) > 0:
                data_center.append(hData.GetBinCenter(i+1))
                data_content.append(hData.GetBinContent(i+1))

        data_center = np.array(data_center)
        data_content = np.array(data_content)


        pdf_center, pdf_content = [], []
        hPdf = TPdf[iPDF]
        for i in range(hPdf.GetNbinsX()):
            pdf_center.append(hPdf.GetBinCenter(i+1))
            pdf_content.append(hPdf.GetBinContent(i+1) * bin_width)

        pdf_center = np.array(pdf_center)
        pdf_content = np.array(pdf_content)
        


        plt.figure(iMass)

        plt.errorbar(data_center, data_content, yerr=np.sqrt(data_content), fmt="o", ms=4, color="blue")
        plt.plot(pdf_center, pdf_content, "-", color='red')
        plt.xlim(-0.02, 0.02)

        iMass+=1

    plt.xlabel("time/s")
    plt.ylabel("#. of Event")
    plt.show()


    """

    bin_width = 0.01
    for iPDF in range(0, nPDF, 1):

        pdf_center, pdf_content = [], []
        hPdf = TPdf[iPDF]
        for i in range(hPdf.GetNbinsX()):
            pdf_center.append(hPdf.GetBinCenter(i+1))
            pdf_content.append(hPdf.GetBinContent(i+1) * bin_width)

        pdf_center = np.array(pdf_center)
        pdf_content = np.array(pdf_content)
        

        plt.plot(pdf_center, pdf_content, "-", label="nuMass%.1feV"%(iPDF/10.))
        plt.xlim(-0.001, 0.01)

    plt.legend()
    plt.show()

    """







