#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: Liangjian Wen ---> wenlj@ihep.ac.cn
# Created Time : Sat Oct 31
# Generate different (time, Ev) PDFs

# Updated by Miao Yu ---> miaoyu@ihep.ac.cn
# Created Time : Mon Aug 9, 2021
"""

import numpy as np
import ROOT
from ROOT import TH1D, TH2D, TCanvas
import scipy.stats as st
from scipy.stats import rv_continuous
from array import array
import sys
import os
import namespace as ns
libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)

if __name__ == "__main__":

    if len(sys.argv) < 10:
        print("Usage: %s [modelNum] [channel] [nuType] [nuMass] [dist] [MH] [Evismin] [Evismax] [no]"  % (sys.argv[0]))
        print("modelNum: check the SNsim/simulation/data directory")
        print("channel : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("nuMass  : neutrino mass in eV unit")
        print("dist    : SN distance in kpc unit")
        print("MH      : mass hierarchy for osc")
        print("Evismin : low edge of Evis" )
        print("Evismax : low edge of Evis" )
        print("No      : current copy no")
        print("python generatePDFs.py 82503 2 1 0.0 1.0 1 0.1 1 5")
        sys.exit(1)
    

    # ------------- Configure Suprenova model ---------------#
    snDet = ROOT.SNdetect.instance()
    # set SN model
    #modelNum = 82503
    modelNum = int( sys.argv[1] )
    snDet.setSrcModel(modelNum)
    modelSrc = snDet.getPointerSrc()
    
    # enum channelName {NuP, NuE, IBD, NCC, BCC, CNC};
    chaName = ['nup','nue','IBD']
    #chaname = snDet.IBD
    chaname = int( sys.argv[2] )
    print("channelName: ", chaname, chaName[chaname])
    snDet.initChannel(chaname)

    dist = float( sys.argv[5] )   # unit: kpc
    print('dist: ', dist, ' kpc')
    modelSrc.setSNDistance(dist)

    Ethr = 0.1
    snDet.getPointerEffectLS().setThresholdE(Ethr)

    snDet.initFCN() # Must call this function, otherwise fEvmax is not set

    nuType = int( sys.argv[3] )
    # -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x
    nuTypeName = ['nu_e', 'anti_nue', 'nu_x']
    print("neutrino type: ", nuType, nuTypeName[nuType])
    
    nuMass = float( sys.argv[4] )# unit: eV
    print('nuMass: ', nuMass, ' eV')

    MH = int( sys.argv[6])  
    print('MH :', MH)
    MHName = ['NoOsc', 'NH', 'IH']

    no = int(sys.argv[9])

    # different mass hierarchy of neutrinos, 0 no osc; 1 NH; 2 IH

    # ------------- Define Histograms ---------------#
    # time binning 
    tmin, tmax = -0.1, 0.1
    t1_thr, t2_thr = -0.025, 0.1
    step_t0, step_t1, step_t2 = 0.00005, 0.00005, 0.0005
    #step_t0, step_t1, step_t2 = 0.001, 0.001, 0.001
    npt0 = round( (t1_thr-tmin)/step_t0 )
    npt1 = round( (t2_thr-t1_thr)/step_t1 )
    npt2 = round( (tmax-t2_thr)/step_t2 )
    nbin_time = npt0 + npt1 + npt2
    
    # time binning
    binning_time = array('d', [])
    for ip in range(npt0):
        binning_time.append( tmin+step_t0*ip )
    
    for ip in range(npt0, npt0 + npt1):
        binning_time.append( t1_thr+step_t1*(ip-npt0) )

    for ip in range(npt0 + npt1, npt0 + npt1 + npt2):
        binning_time.append( t2_thr+step_t2*(ip-npt0-npt1) )
    
    binning_time.append(tmax)
    # for debug
    print('npt0, start, end: ', npt0, binning_time[0], binning_time[npt0])
    print('npt1, start, end: ', npt1, binning_time[npt0], binning_time[npt0+npt1])
    print('npt2, start, end: ', npt2, binning_time[npt0+npt1], binning_time[npt0+npt1+npt2])
    print('len(binning_time), nbin_time, last bin:', len(binning_time), nbin_time, binning_time[npt0+npt1+npt2])

    # energy binning
    Evismin, Evismax, step_Evis = Ethr, 1.0, 0.01
    Evmin,   Evmax,   step_Ev   = 0.0, 80.0, 0.5
    
    if chaname == snDet.NuP:
        Evismax, step_Evis = 5.0, 0.01
    if chaname == snDet.NuE:
        Evismax, step_Evis = 25,  0.05
        #Evismax, step_Evis = 25,  0.5
        Evmax,   step_Ev   = 60, 0.2
        Evismin = float( sys.argv[7] )
        Evismax = float( sys.argv[8] )
        print("Evis range [ %.2f , %.2f ]" %(Evismin, Evismax))
    if chaname == snDet.IBD:
        Evismin, Evismax, step_Evis = 1.02, 25+1.02,  0.02
        Evmin, Evmax, step_Ev = 1.02+1.3-0.51, 30+1.02+1.3-0.51, 0.02

    nbins_Ev   = round( (Evmax-Evmin)/step_Ev )
    nbins_Evis = round( (Evismax-Evismin)/step_Evis )
    print('Evismin, Evismax, nbins_Evis', Evismin, Evismax, nbins_Evis)
    print('Evmin, Evmax, nbins_Ev', Evmin, Evmax, nbins_Ev)


    # ------------- Output PDFs ---------------#
    path = "./etSpec/fineSpec/"
    imod = array('i', [modelNum])
    icha = array('i', [chaname])
    tStart, tEnd = array('d', [-999.]), array('d', [-999.])
    outFileName = ns.pdfFileName(modelNum, dist, chaname, nuMass, MH, nuType, Evismax, no)
    print(outFileName)
    fout = ROOT.TFile(outFileName, "recreate")

    # ------------- TTree Info ---------------#
    tinfo = ROOT.TTree("tinfo", "model information")
    tinfo.Branch("model", imod, "model/I")
    tinfo.Branch("channel", icha, "channel/I")
    tinfo.Branch("tStart", tStart, "tStart/D")
    tinfo.Branch("tEnd", tEnd, "tEnd/D")
    modelSrc.getTimeRange(tStart, tEnd, chaname)
    tinfo.Fill()


    # ------------- 2D distribution ---------------#
    h2d_etSpec, h2d_eNuTSpec, h2d_etSpecBKG = [], [], []

    evisBins = array('d', [0.2, 0.4, 0.6, 1.0, 4.0, 6.0, 8.0, 10.0])
    dummyEv   = ROOT.TH1D('dummyEv', 'dummyEv', nbins_Ev, Evmin, Evmax)
    dummyEvis = ROOT.TH1D('dummyEvis', 'dummyEvis', len(evisBins)-1, evisBins)
    dummyT    = ROOT.TH1D('dummyT', 'dummyT', nbin_time, binning_time)

    histname = ns.pdfEvisHist2DName(modelNum, chaname, MH)
    histTitle = histname
    histTmp = ROOT.TH2D(histname, histTitle, nbin_time, binning_time, nbins_Evis, Evismin, Evismax)
    h2d_etSpec.append( histTmp )
    h2d_etSpec[0].GetXaxis().SetTitle("time [s]")
    h2d_etSpec[0].GetYaxis().SetTitle("E_{vis} [MeV]")
    h2d_etSpec[0].GetZaxis().SetTitle('%s event'%chaName[chaname])

    histname = ns.pdfBkgHist2DName(modelNum, chaname, MH)
    histTitle = histname
    histTmp = ROOT.TH2D(histname, histTitle, nbin_time, binning_time, nbins_Evis, Evismin, Evismax)
    h2d_etSpecBKG.append( histTmp )
    h2d_etSpecBKG[0].GetXaxis().SetTitle("time [s]")
    h2d_etSpecBKG[0].GetYaxis().SetTitle("E_{vis} [MeV]")
    h2d_etSpecBKG[0].GetZaxis().SetTitle('%s event'%chaName[chaname])

    histname = ns.pdfEvHist2DName(modelNum, chaname, MH)
    histTmp = ROOT.TH2D(histname, histTitle, nbin_time, binning_time, nbins_Ev, Evmin, Evmax)
    h2d_eNuTSpec.append( histTmp )
    h2d_eNuTSpec[0].GetXaxis().SetTitle("time [s]")
    h2d_eNuTSpec[0].GetYaxis().SetTitle("E_{v} [MeV]")
    h2d_eNuTSpec[0].GetZaxis().SetTitle('%s event'%chaName[chaname])

    
    
    # ----------- ROI ------------- #
    tBurstEnd = 1.0
    binBurstEnd = h2d_etSpec[0].GetXaxis().FindBin(tBurstEnd)
    if binBurstEnd>h2d_etSpec[0].GetNbinsX():
        binBurstEnd = h2d_etSpec[0].GetNbinsX()
    print('tBurstEnd, binBurstEnd: ', h2d_etSpec[0].GetXaxis().GetBinCenter(binBurstEnd), binBurstEnd)


    # ------------- Modify the PDFs due to different neutrino masses ---------------#
    peffLS = snDet.getPointerEffectLS()
    Npar = 0
    # number of protons
    if chaname == snDet.NuP or chaname == snDet.IBD:
        Npar = peffLS.getNumberOfProton()
    # number of electrons
    if chaname == snDet.NuE:
        Npar = peffLS.getNumberOfProton()+6*peffLS.getNumberOfCarbon()
    # number of carbon-12, 98.9% of carbons
    if chaname == snDet.NCC or chaname == snDet.BCC or chaname == snDet.CNC:
        Npar = 0.989*peffLS.getNumberOfCarbon()
    print('Npar', Npar)

    timeTmp = array('d', [-999.]) # unit: s
    #me = 0.511; #MeV

    nCalc = 0
    Nabnormal = 0
    for ipt in range(binBurstEnd):
        itwidth = binning_time[ipt+1] - binning_time[ipt]

        for k in range(nbins_Evis):
            EvisTmp = Evismin + (k+0.5) * step_Evis

            if nCalc%1000 == 0:
                print("nCalc : ", nCalc)
            nCalc += 1


            if chaname == snDet.NuE:
                for j in range(nbins_Ev):
                    EvTmp = Evmin + (j+0.5)*step_Ev
                    timeTmp[0] = binning_time[ipt] + 0.5*itwidth
                    
                    # single flavour
                    if nuType != -1 :
                        fluence = modelSrc.snFluenceDetAtTime(timeTmp, nuMass, EvTmp, nuType, MH)
                        dxs = peffLS.differentialXS(EvTmp, EvisTmp, nuType)
                        if fluence > 0. and dxs >0.:

                            h2d_etSpec[0].Fill(timeTmp[0], EvisTmp, Npar*fluence*dxs*step_Ev)

                            h2d_eNuTSpec[0].Fill(timeTmp[0], EvTmp, Npar*fluence*dxs)


                    if nuType == -1 :
                        tot_fluence = 0
                        for i in range(6):
                            timeTmp[0] = binning_time[ipt] + 0.5*itwidth    # set time repeatedly due to each time the time will be updated during snFluenceDetAtTime
                            fluence = modelSrc.snFluenceDetAtTime(timeTmp, nuMass, EvTmp, i, MH)
                            dxs = peffLS.differentialXS(EvTmp, EvisTmp, i)
                            tot_fluence += fluence*Npar*dxs
                            #print(timeTmp[0], EvTmp, EvisTmp, fluence, dxs, tot_fluence)
                            
                        if tot_fluence > 0. :
                            h2d_etSpec[0].Fill(timeTmp[0], EvisTmp, tot_fluence*step_Ev)
                            h2d_eNuTSpec[0].Fill(timeTmp[0], EvTmp, tot_fluence)



    fout.cd()
    tinfo.Write()

    # ----------------- Background ------------------ #
    tmpBKG = 60/86400./15
    # Set background from reactors and others
    Xbins = h2d_etSpecBKG[0].GetNbinsX()
    Ybins = h2d_etSpecBKG[0].GetNbinsY()
    for i in range(1, Xbins+1):
        Xwidth = h2d_etSpecBKG[0].GetXaxis().GetBinWidth(i)

        for j in range(1, Ybins+1):
            if h2d_etSpecBKG[0].GetYaxis().GetBinCenter(j)<1.0:
                continue

            Ywidth = h2d_etSpecBKG[0].GetYaxis().GetBinWidth(j)
            #val = h2d_etSpec[0].GetBinContent(i, j)

            h2d_etSpecBKG[0].SetBinContent(i, j, tmpBKG)
    

    h2d_etSpec[0].Write()
    h2d_eNuTSpec[0].Write()
    h2d_etSpecBKG[0].Write()



    fout.Close()
