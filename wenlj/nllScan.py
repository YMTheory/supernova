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
    return (x1, x2, x3)

def func(x, a0, t0, a1):
    return a0*np.power(x-t0, 2) + a1

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
    
    nuMass1 = 0.0
    #nuMass2 = 0.2
    modelNum = 82503
    chaname = 1
    chaName = ['nup','nue','IBD']
    dist = 3.0
    MH = 0
    mNu, minNLL = [], []

    prefix = "./etSpec/TEvisPDF_"
    prefix += 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)
    
    with PdfPages('nllProfile_data%2.1feV_%2.1fkpc_1DFit.pdf'%(nuMass1, dist)) as pdf:
        for k in range(11, 12):
            nuMass2 = k * 0.1
            fname = prefix+'_data%2.1feV_pdf%2.1feV_%2.1fkpc_1DFitRaw.txt'%(nuMass1, nuMass2, dist)
            if os.path.exists(fname) == False:
                continue
            
            deltaT, nllVal, locMinNLL = getNLL(fname)
            #popt, pcov = curve_fit(func, deltaT, nllVal)
            minNLL.append(locMinNLL)
            mNu.append(nuMass2)
                
            plotx = np.arange(-0.2, 0.2)
            plt.plot(deltaT, nllVal,  "o", ms=5, label="Sinnock 1969")
            #plt.plot(plotx, func(plotx, *popt), 'r-',
            #   label='fit: a=%.3f, t0=%.4f, minNLL=%.6f' % tuple(popt))
            plt.xlabel('$\delta_{T}$ (s)')
            plt.ylabel('NLL')
            titlename = 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)
            titlename += '_data%2.1feV_pdf%2.1feV'%(nuMass1, nuMass2)
            plt.title(titlename)
            pdf.savefig()
            plt.close()

        minNLL = np.array(minNLL)
        mNu = np.array(mNu)
        print(minNLL)
        plotx = np.arange(-0.2, 0.2)
        plt.plot(mNu, minNLL,  "o", ms=5, label="Sinnock 1969")
        plt.xlabel(r'$m_{\nu}$ (eV)')
        plt.ylabel('local minNLL')
        titlename = 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)
        titlename += '_data%2.1feV'%(nuMass1)
        plt.title(titlename)
        pdf.savefig()
        plt.close()
