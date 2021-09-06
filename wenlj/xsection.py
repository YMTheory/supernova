#!/usr/bin/env python
# -*- coding=utf8 -*-

import numpy as np
import ROOT
import sys
import os
import matplotlib.pyplot as plt
libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)

if __name__ == "__main__":

    # ------------- Configure Suprenova model ---------------#
    snDet = ROOT.SNdetect.instance()
    # set SN model
    modelNum = 82503
    snDet.setSrcModel(modelNum)
    modelSrc = snDet.getPointerSrc()
    
    # enum channelName {NuP, NuE, IBD, NCC, BCC, CNC};
    chaName = ['nup','nue','IBD']
    chaname = 1
    snDet.initChannel(chaname)

    dist = 10.0
    modelSrc.setSNDistance(dist)

    Ethr = 0.1
    snDet.getPointerEffectLS().setThresholdE(Ethr)

    snDet.initFCN() # Must call this function, otherwise fEvmax is not set

    peffLS = snDet.getPointerEffectLS()

    totXS_eES_nue, totXS_eES_nuebar, totXS_eES_nux, totXS_eES_nuxbar = [], [], [], []
    Enu = np.arange(0.0, 60.0, 0.1) 
    for i in Enu:
        dxs = peffLS.totalXS(i, 0)
        totXS_eES_nue.append(dxs)
        dxs = peffLS.totalXS(i, 1)
        totXS_eES_nuebar.append(dxs)
        dxs = peffLS.totalXS(i, 2)
        totXS_eES_nux.append(dxs)
        dxs = peffLS.totalXS(i, 3)
        totXS_eES_nuxbar.append(dxs)

    chaname = 0
    snDet.initChannel(chaname)
    snDet.initFCN() # Must call this function, otherwise fEvmax is not set
    peffLS = snDet.getPointerEffectLS()
    totXS_pES = []
    for i in Enu:
        dxs = peffLS.totalXS(i, 0)
        totXS_pES.append(dxs)

    plt.plot(Enu, totXS_eES_nue, label=r"$eES \nu_e$")
    plt.plot(Enu, totXS_eES_nuebar, label=r"$eES \bar\nu_e$")
    plt.plot(Enu, totXS_eES_nux, label=r"$eES \nu_x$")
    plt.plot(Enu, totXS_eES_nuxbar, label=r"$eES \bar\nu_x$")

    plt.plot(Enu, totXS_pES, "--", label=r"$pES$")

    plt.semilogy() 
    plt.legend()
    plt.xlabel(r"$E_{\nu}/MeV$")
    plt.ylabel(r"$\sigma/cm^{-2}$")
    plt.savefig("xsection_eES+pES.pdf")
    plt.show()
