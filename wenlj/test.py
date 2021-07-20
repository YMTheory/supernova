#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Thu Aug 27 20:22:38 2020
# File Name: hist.py
# Example for histograms in matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gSystem.Load("/mnt/c/Users/LiangjianWen/Documents/JUNO/Physics/SNsim/simulation/lib/libSNsim.so")
from ROOT import TH1D, TH2D, TCanvas
import scipy.stats as st
from scipy.stats import rv_continuous
from array import array

#plt.style.use('fivethirtyeight')
#plt.style.use('bmh')


def draw_hist1d():
    """ draw a 1d histogram """
    data = np.random.normal(100, 20, 10000)
    plt.hist(data, bins=100, range=(0, 200), histtype="step", color="royalblue")

def draw_hist2d():
    datax = np.random.normal(100, 20, 10000)
    datay = np.random.normal(100, 30, 10000)
    ec = plt.hist2d(datax, datay, bins=(100, 100), range=((0, 200), (0, 200)))
    cbar = plt.colorbar()
    cbar.set_label('entries')

class nuT_gen(rv_continuous):
    "generate the time of radiated neutrinos"
    def _pdf(self, x, *args):
        modelSrc = args[0]
        nuType = args[1]
        MH = args[2]
        return modelSrc.oneSNFluenceDetTimeIntegE(x, nuType, MH)

#def fcnEv2Time(self, type, MH):
    # Convert Ev to visible energy for different channels
#    xs = snDet.getPinterEffectLS().differentialXS(double E, double T, int type)

if __name__ == "__main__" :
#with plt.style.context("Style/Paper.mplstyle"):

#    fig = plt.figure(figsize=(10, 4))
    # PLOT 1
#    ax = fig.add_subplot(1, 2, 1)
#    draw_hist1d()
    # PLOT 2
#    ax = fig.add_subplot(1, 2, 2)
#    draw_hist2d()
    
#    plt.show()
    snDet = ROOT.SNdetect.instance()
    # set SN model
    snDet.setSrcModel(82503)
    # enum channelName {NuP, NuE, IBD, NCC, BCC, CNC};
    snDet.initChannel(snDet.NuE)
    print("channelName", snDet.NuE)
    snDet.initFCN() # Must call this function, otherwise fEvmax is not set

    nuType = 0   # -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x
    nuTime = 0.2  # unit: s
    nuEvis = 10.1 # unit: MeV
    MH = 0       # different mass hierarchy of neutrinos, 0 no osc; 1 NH; 2 IH
    nuNum = []
    snDist, tMin, tMax = 10, array('d', [0.]), array('d', [0.])



    #snDet.fcnAnaEv2TatTime(1, 0, 0.5, nuType, MH)

    snModel = snDet.getPointerSrc()
    nuTGen = nuT_gen(name="nuTGen")
    snModel.getTimeRange(tMin, tMax, nuType)
    print ("tMin, tMax: ", tMin[0], tMax[0])
    st.rvs_ratio_uniforms(nuTGen.pdf(snModel, nuType, MH), 0.5, tMin[0], tMax[0], size=10)


    for i in range(5):
        nuTime = i*0.005

        # get neutrino spectrum at given time
    #    snDet.getPointerSrc().oneSNFluenceDetAtTime(nuTime, nuE, nuType)

        # Convert Ev to visible energy for different channels


        # add time delay of arrvial neutrinos due to non-zero neutrino mass
    #    DeltaT = 5.14 * (mNu*mNu) * (100.0/nuEv/nuEv) * （snDist/10)
        
        # Fill (deltaT, Evis)


    #    if i%10==0:
    #        print(nuTime, snDet.getEvisSpectrumAtTime(nuTime, nuEvis, nuType, MH))
        #nuNum.append(snDet.getEvisSpectrumAtTime(nuTime, nuEvis, nuType, MH))
    
    print(nuNum)

    #fig2 = plt.figure(figsize=(5,5))
    #ax = fig2.add_subplot()
    ##ec = plt.hist(relTobsNuebar, 50, histtype='stepfilled', facecolor='g', range=(-1,9))
    #ec = plt.hist(nuNum, 50, histtype='step', facecolor='b', range=(0,5))
    #ax.set_xlabel('E_{vis} (MeV)')
    #plt.show()
    #print("test")
