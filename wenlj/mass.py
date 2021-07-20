#!/usr/bin/env python

import matplotlib.pyplot as plt
import ROOT
import numpy as np

# plot with function or data
def PlotAFCN(nPts=200):
    x=np.linspace(0.1,10,nPts)
    y=np.sin(x)/x
    plt.plot(x,y,label='sin(x)/x',color='coral',linestyle='-.')
    plt.legend(loc='best')
    plt.show()

def calcDeltaT(m_nu, e_nu, D): # unit: ms
    DeltaT = 5.14 * m_nu*m_nu / (e_nu*e_nu) * (D/10)
    return DeltaT 

def CalcEvis(opt, e_nu):
    if opt == 'eES':
        Evis = e_nu
    elif opt == 'IBD':
        Evis = e_nu - 0.8
    elif opt == 'pES':
        Evis = e_nu
    else:
        Evis = -1.

    return Evis

if __name__ == "__main__":
    #Global plot style option
    # plt.style.use("Style/Paper.mplstyle")
    
    #Plot style within given scope
    #with plt.style.context("Style/Paper.mplstyle"):
    #    PlotAFCN()
        #PlotByRoot()
        #PlotWithRootData()

    ROOT.gSystem.Load("/junofs/users/wenlj/juno/supernova/SNsim/simulation/lib/libSNsim.so")
    o = ROOT.SNdetect.instance()

    mNu = 0.2 # unit: 1 eV
    dist = 10 # unit: kpc
    Nevt = 10000

    eNu = np.random.normal(2, 0.2, Nevt)
    #Tgen = np.random.normal(10, 2, Nevt) # unit: ms
    Tgen = np.random.triangular(0.0, 5.0, 5.0, Nevt) # unit: ms
    EvisNuebar = []
    TobsNuebar = []
    relTobsNuebar = []
    relTgenNuebar = []

    #print (eNu)

    # generate (E, t) from SN spectra
    for i in range(Nevt):
        #print (eNu[i])

        dT = calcDeltaT(mNu, eNu[i], dist)
        tObs = Tgen[i] + dT

        Evis = CalcEvis('IBD', eNu[i])
        EvisNuebar.append(Evis)

        TobsNuebar.append(tObs)

    tObsMin = np.min(TobsNuebar)
    tGenMin = np.min(Tgen)

    # subtract the earliest neutrino observation time
    for i in range(Nevt):
        relTobsNuebar.append( TobsNuebar[i] - tObsMin)
        relTgenNuebar.append( Tgen[i] - tGenMin)

    #np.set_printoptions(formatter={'float':'{:0.3f}', format}}
    #print (TobsNuebar)
    #print (Tgen)
    #print (tObsMin, tGenMin)
    #print (EvisNuebar)
    #print (TobsNuebar)

    #fig, ax = plt.subplots(figsize=(6,5))
    #ec = plt.hist2d(EvisNuebar, TobsNuebar, bins=(100, 100), range=((0, 10), (0, 10)))
    #cbar = plt.colorbar()
    #ax.set_xlabel('Tobs (ms)')
    #ax.set_ylabel('Evis (MeV)')

    fig = plt.figure(figsize=(11,5))
    ax = fig.add_subplot(1, 2, 1)
    #ec = plt.hist2d(TobsNuebar, EvisNuebar, bins=(100, 100), range=((0, 10), (0, 10)))
    ec = plt.hist2d(relTobsNuebar, EvisNuebar, bins=(100, 100), range=((-1, 9), (0, 4)))
    cbar = plt.colorbar()
    ax.set_xlabel('Tobs (ms)')
    ax.set_ylabel('Evis (MeV)')

    ax = fig.add_subplot(1, 2, 2)
    #ec = plt.hist2d(Tgen, EvisNuebar, bins=(100, 100), range=((0, 10), (0, 10)))
    ec = plt.hist2d(relTgenNuebar, EvisNuebar, bins=(100, 100), range=((-1, 9), (0, 4)))
    cbar = plt.colorbar()
    ax.set_xlabel('Tgen (ms)')
    ax.set_ylabel('Evis (MeV)')

    fig2 = plt.figure(figsize=(5,5))
    #ec = plt.hist(relTobsNuebar, , density=True, facecolor='b', alpha=0.75)
    ax = fig2.add_subplot()
    ec = plt.hist(relTobsNuebar, 50, histtype='stepfilled', facecolor='g', range=(-1,9))
    ec = plt.hist(relTgenNuebar, 50, histtype='step', facecolor='b', range=(-1,9))
    ax.set_xlabel('T (ms)')

    plt.show()

