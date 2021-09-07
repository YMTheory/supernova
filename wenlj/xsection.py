#!/usr/bin/env python
# -*- coding=utf8 -*-

import numpy as np
import ROOT
import sys
import os
import matplotlib.pyplot as plt
libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)

def pES_XS(Ev, Tp, nutype):  # all in MeV
    Mp = 938.    # MeV
    index = 7.88e-42   # cm2/MeV
    cV, cA = 0.04, 0.635
    if nutype == 0 :
        return index * ( (cV+cA)**2*Ev**2 + (cV-cA)**2*(Ev-Tp)**2 - (cV**2-cA**2)*Mp*Tp) / Ev**2
    if nutype == 1:
        return index * ( (cV-cA)**2*Ev**2 + (cV+cA)**2*(Ev-Tp)**2 - (cV**2-cA**2)*Mp*Tp) / Ev**2



def generator():

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


def TpMax(Ev):
    Mp = 938
    return 2*Ev**2/(Mp+2*Ev)


def calculation():

    Env = 30
    Tmax = TpMax(Env)
    print("%.1f MeV neutrino -> maximum proton KE = %.2f MeV" %(Env, Tmax))

    Tp = np.arange(0, Tmax, 0.1)

    xs_nu, xs_antinu = [], []
    for i in Tp:
        xs_nu.append(pES_XS(Env, i, 0))
        xs_antinu.append(pES_XS(Env, i, 1))



    plt.plot(Tp, xs_nu, "-", color="blue", label="Enu=30MeV, nu")
    plt.plot(Tp, xs_antinu, "--", color="blue", label="Enu=30MeV, anti-nu")

    Env = 40
    Tmax = TpMax(Env)
    print("%.1f MeV neutrino -> maximum proton KE = %.2f MeV" %(Env, Tmax))

    Tp = np.arange(0, Tmax, 0.1)

    xs_nu, xs_antinu = [], []
    for i in Tp:
        xs_nu.append(pES_XS(Env, i, 0))
        xs_antinu.append(pES_XS(Env, i, 1))



    plt.plot(Tp, xs_nu, "-", color="seagreen", label="Enu=40MeV, nu")
    plt.plot(Tp, xs_antinu, "--", color="seagreen", label="Enu=40MeV, anti-nu")
    Env = 50
    Tmax = TpMax(Env)
    print("%.1f MeV neutrino -> maximum proton KE = %.2f MeV" %(Env, Tmax))

    Tp = np.arange(0, Tmax, 0.1)

    xs_nu, xs_antinu = [], []
    for i in Tp:
        xs_nu.append(pES_XS(Env, i, 0))
        xs_antinu.append(pES_XS(Env, i, 1))



    plt.plot(Tp, xs_nu, "-", color="orange", label="Enu=50MeV, nu")
    plt.plot(Tp, xs_antinu, "--", color="orange", label="Enu=50MeV, anti-nu")

    plt.legend()
    plt.xlabel(r"$T_p/MeV$")
    plt.ylabel(r"$\frac{d\sigma}{dT_p} [cm^{2}/MeV]$")

    plt.savefig("pES_xsection.pdf")
    plt.show()





if __name__ == "__main__" :
    calculation()
