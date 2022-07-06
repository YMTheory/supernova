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
    
    mod = 82503
    MH = 1

    if len(sys.argv) == 3:
        mod = int(sys.argv[1])
        MH  = int(sys.argv[2])

    # declare variables
    Tmin, Tmax = 0.003, 0.08     # Fitting time range
    time   = ROOT.RooRealVar("time", "time", Tmin, Tmax)
    time.setRange('nllRange', Tmin, Tmax)
    

    # load dataset :
    cha = 0
    filename = "/junofs/users/miaoyu/supernova/production/Data/10kpc/data_mod%d_cha1_mh%d.root"%(mod, MH)
    outpdfname = "./outputs/fit_mod%d_cha%_mh%d.pdf" %(mod, cha, MH)

    print(filename)
    dataFile  = ROOT.TFile(filename, 'read')
    inputTree = dataFile.Get('evt')
    evtID    = ROOT.RooRealVar("evtID", "evtID", 0, 500)
    newArgSet = ROOT.RooArgSet(time)

    
    nStat = 500




