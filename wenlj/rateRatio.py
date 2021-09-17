import ROOT
import matplotlib.pyplot as plt
import sys
import numpy as np
import namespace as ns


def getIntegral(hist2d, x1, x2, y1, y2, opt):
    binx1 = hist2d.GetXaxis().FindBin(x1)
    binx2 = hist2d.GetXaxis().FindBin(x2)
    biny1 = hist2d.GetYaxis().FindBin(y1)
    biny2 = hist2d.GetYaxis().FindBin(y2)

    return hist2d.Integral(binx1, binx2, biny1, biny2, opt)


if __name__ == "__main__"  :

    if len(sys.argv) < 8:
        print("Usage: %s [modelNum] [channel] [nuType] [nuMass] [pdfMH] [dist] [Ethr] [group]  "%(sys.argv[0]))
        print("modelNum: check the SNsim/simulation/data directory")
        print("channel : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("nuMass  : neutrino mass that used to generate dataset, in eV unit")
        print("dist    : SN distance in kpc unit")
        print("Ethr    : fit Energy threshold")
        print("group   : index of the groups")
        print("example : python nllFit.py 82503 1 0 0.0 2 0.2 1 10.0 0.1 1")
        sys.exit(1)

    ##################################
    # Input variables
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
    print('nuMass for pdf : ', nuMass, ' eV')

    dist = float( sys.argv[5] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[6] )# unit: kpc
    print('fit Ethr: ', Ethr, ' MeV')

    group = int( sys.argv[7] )# group number
    print('group number: ', group)

    # ###################################
    # # Obtain the fitting PDF
    filename = ns.pdfSumFileName(modelNum, dist, chaname, nuMass, 1, nuType, 25.00)
    print(filename)
    pdfFileNO  = ROOT.TFile(filename, 'read')
    inputHistNO = pdfFileNO.Get("TEvisPdf")
    #inputBKG  = pdfFile.Get("TBkgEvisPdf")
    filename = ns.pdfSumFileName(modelNum, dist, chaname, nuMass, 2, nuType, 25.00)
    print(filename)
    pdfFileIO  = ROOT.TFile(filename, 'read')
    inputHistIO = pdfFileIO.Get("TEvisPdf")


    hist1 = ROOT.TH2D("NO", "", 10, 0.015, 0.02, 40, 0.04, 0.06)
    hist2 = ROOT.TH2D("IO", "", 10, 0.015, 0.02, 40, 0.04, 0.06)
    
    Tmin = -0.015
    i = 0
    for x in np.arange(0.015, 0.02, 0.0005) :
        i += 1
        j = 0
        for y in np.arange(0.04, 0.06, 0.0005) :
            j += 1
            num1 = getIntegral(inputHistNO, Tmin, x, 0.1, 10, "width")
            num2 = getIntegral(inputHistNO, Tmin, y, 0.1, 10, "width")
            hist1.SetBinContent(i, j, (num2-num1)/num2)
            num1 = getIntegral(inputHistIO, Tmin, x, 0.1, 10, "width")
            num2 = getIntegral(inputHistIO, Tmin, y, 0.1, 10, "width")
            hist2.SetBinContent(i, j, (num2-num1)/num2)


    out = ROOT.TFile("out.root", "recreate")
    hist1.Write()
    hist2.Write()
    out.Close()



