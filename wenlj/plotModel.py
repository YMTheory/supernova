import ROOT
import math
import string
import matplotlib.pyplot as plt
import sys, os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

def getValue(strline):
    values = strline.split()
    x1 = float(values[0])
    x2 = float(values[1])
    x3 = float(values[2])
    x4 = float(values[3])
    return (x1, x2, x3, x4)

def func(x, a0, t0, a1):
    return a0*np.power(x-t0, 2) + a1

def getData(filename):
    nuTime, lum, nuE, nuEsquare = [], [], [], []
    
    if os.path.exists(filename) == False:
        return False

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    print('number of lines', len(lines))

    for k in range(len(lines)):
        (x1, x2, x3, x4) = getValue( lines[k] )
        nuTime.append(x1)
        lum.append(x2)
        nuE.append(x3)
        nuEsquare.append(x4)
    return nuTime, lum, nuE, nuEsquare
    # deltaT = np.array(deltaT)
    # nllVal = np.array(nllVal)
    # print(deltaT)
    # print(nllVal)
    # locMinNLL = np.min(nllVal)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: %s [modelType] [modelNum] [channel] [nuType]"  % (sys.argv[0]))
        print("modelType: check the SNsim/simulation/data directory")
        print("modelNum: check the SNsim/simulation/data directory")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("example: python plotModel.py 82503 1 0")
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

    nuType = int( sys.argv[2] )
    # -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x
    nuTypeName = ['nu_e', 'nubar_e', 'nu_x']
    print("neutrino type: ", nuType, nuTypeName[nuType])

    modelList = []
    if modelType == 'Garching':
        modelList = [71200,
                 81120,81761,82504,84004,91760,92502,94003,
                 71380,81121,81780,82700,91120,91761,92503,
                 71780,81122,82000,82701,91121,91780,92700,
                 72000,81123,82060,82702,91122,91781,92701,
                 72060,81124,82061,82703,91123,92000,92702,
                 72500,81200,82500,84000,91200,92060,92703,
                 73500,81500,82501,84001,91201,92061,94000,
                 73600,81501,82502,84002,91500,92500,94001,
                 74000,81760,82503,84003,91501,92501,94002 ]

    elif modelType == 'Japan':
        print('process Japan model')
        sys.exit(1)

    pdfFilename = './spectra/model'
    if modelNum != -1:
        modelList = [modelNum]
        pdfFilename = pdfFilename + '%d.pdf'%(modelNum)
    else:
        pdfFilename = pdfFilename + '_all.pdf'

    with PdfPages(pdfFilename) as pdf:
        for num in modelList:
            print(num)
            filename = "./simulation/data/Garching/%d/timedata/neutrino_signal_%s"%(num, nuTypeName[nuType])

            nuTime, lum, nuE, nuEsquare = getData(filename)

            fig, axs = plt.subplots(2, 1, figsize=(12, 8))
            plt.subplot(211)
            plt.plot(nuTime, nuE, '.')
            plt.axis([-0.1,0.1,2.,18.])
            #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
            plt.ylabel('Energy (MeV)')
            plt.xlabel('Time (s)')
            titlename = 'model %d'%(num)
            plt.title(titlename)

            plt.subplot(212)
            plt.plot(nuTime, lum, '.')
            plt.axis([-0.1,0.1,0.,400.])
            #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
            plt.ylabel('Luminosity')
            plt.xlabel('Time (s)')

            pdf.savefig(dpi=400,bbox_inches='tight')
            #plot.savefig('./spectra/model%d.png'%(modelNum),dpi=400,bbox_inches='tight')
            plt.close()
