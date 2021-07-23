#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: Liangjian Wen ---> wenlj@ihep.ac.cn
# Created Time : Sat Oct 31
# Generate different (time, Ev) PDFs
# Modified by Miao Yu ---> miaoyu@ihep.ac.cn for reproduction
"""

import numpy as np
#import matplotlib.pyplot as plt
import ROOT
from ROOT import TH1D, TH2D, TCanvas
import scipy.stats as st
from scipy.stats import rv_continuous
from array import array
import sys
import os
import SNnueXS as nuexs
import SNGarchingIntegFcn as gar
from SNnumGarchingSrc import SNnumGarchingSrc
from detector import detector

if __name__ == "__main__":

    if len(sys.argv) < 7:
        print("Usage: %s [modelNum] [channel] [nuType] [nuMass] [dist] [MHtype]"  % (sys.argv[0]))
        print("modelNum: check the SNsim/simulation/data directory")
        print("channel : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("nuMass  : neutrino mass in eV unit")
        print("dist    : SN distance in kpc unit")
        print("python generatePDFs.py 82503 2 1 0.0 1.0")
        sys.exit(1)

    # ------------- Configure Suprenova model ---------------#
    # set SN model
    print(" ==============> Configurations of Suprenova <=============")
    modelNum = int( sys.argv[1] )
    print("---> Garching model %d" %modelNum)
    
    # enum channelName {NuP, NuE, IBD, NCC, BCC, CNC};
    chaName = ['nup','nue','IBD']
    chaname = int( sys.argv[2] )
    print("---> channelName: ", chaname, chaName[chaname])

    dist = float( sys.argv[5] )   # unit: kpc
    print('---> dist: ', dist, ' kpc')

    Ethr = 0.1

    nuType = int( sys.argv[3] )
    # -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x
    nuTypeName = ['nu_e', 'anti_nue', 'nu_x']
    print("---> neutrino type: ", nuType, nuTypeName[nuType])
    
    nuMass = float( sys.argv[4] )# unit: eV
    print('---> nuMass: ', nuMass, ' eV')
    MHName = ['NoOsc', 'NH', 'IH']
    MH = int( sys.argv[6] )
    # different mass hierarchy of neutrinos, 0 no osc; 1 NH; 2 IH
    print('---> Mass ordering : %s'%MHName[MH])

    SNGar = SNnumGarchingSrc(modelNum, dist)
    det = detector()


    # ------------- Define Histograms ---------------#
    # tmin, tmax = -0.05, 1.0
    # t1_thr, t2_thr = 0.1, 0.6
    # step_t0, step_t1, step_t2 = 0.0005, 0.01, 0.04
    
    tmin, tmax = -0.1, 0.1
    #t1_thr, t2_thr = -0.025, 0.05
    t1_thr, t2_thr = -0.025, 0.1
    step_t0, step_t1, step_t2 = 0.00005, 0.00005, 0.0005
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
    
    if chaname == 0:
        Evismax, step_Evis = 5.0, 0.01
    if chaname == 1:
        Evismax, step_Evis = 25,  0.05
        Evmax,   step_Ev   = 60, 0.2
        # Evismax, step_Evis = 8,  0.5
        # Evmax,   step_Ev   = 80, 1.0
    if chaname == 2:
        #Evismax, step_Evis = 5.0, 0.01
        Evismin, Evismax, step_Evis = 1.02, 25+1.02,  0.02
        Evmin, Evmax, step_Ev = 1.02+1.3-0.51, 30+1.02+1.3-0.51, 0.02

    nbins_Ev   = round( (Evmax-Evmin)/step_Ev )
    nbins_Evis = round( (Evismax-Evismin)/step_Evis )
    print('Evismin, Evismax, nbins_Evis', Evismin, Evismax, nbins_Evis)
    print('Evmin, Evmax, nbins_Ev', Evmin, Evmax, nbins_Ev)

    # ------------- Output PDFs ---------------#
    path = "./outputs"
    imod = array('i', [modelNum])
    icha = array('i', [chaname])
    tStart, tEnd = array('d', [-999.]), array('d', [-999.])
    #fout = ROOT.TFile('%s/TEvisPDF_mod%d_cha%d%s_mh%d_mNu%2.1feV_%2.1fkpc_0.1s_v4.root'
    fout = ROOT.TFile('%s/TEvisPDF_mod%d_cha%d%s_mh%d_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25.root'
            %(path, modelNum, chaname, chaName[chaname], MH, nuMass, dist), 'recreate')
    tinfo = ROOT.TTree("tinfo", "model information")
    tinfo.Branch("model", imod, "model/I")
    tinfo.Branch("channel", icha, "channel/I")
    tinfo.Branch("tStart", tStart, "tStart/D")
    tinfo.Branch("tEnd", tEnd, "tEnd/D")
    tinfo.Fill()


    # ------------- 2D distribution ---------------#
    nMH = 3
    h2d_etSpec, h2d_eNuTSpec, h1d_eNuSpec, h1d_eEvisSpec = [], [], [], []
    h2d_etSpecBKG, h2d_eNuTSpecSlice = [], []
    #h2d_eNuTCount, h2d_etCount = [], []
    
    evisBins = array('d', [0.2, 0.4, 0.6, 1.0, 4.0, 6.0, 8.0, 10.0])
    dummyEv = ROOT.TH1D('dummyEv', 'dummyEv', nbins_Ev, Evmin, Evmax)
    dummyEvis = ROOT.TH1D('dummyEvis', 'dummyEvis', len(evisBins)-1, evisBins)
    dummyT = ROOT.TH1D('dummyT', 'dummyT', nbin_time, binning_time)
    for imh in range(nMH):
        histname = 'hET_mod%d_cha%d_mh%d'%(modelNum, chaname, imh)
        histTitle = '(T_{obs}, E_{vis}) for '
        histTitle += 'model%d, %s, %s, m_{#nu}%2.1f eV'%(modelNum, chaName[chaname], MHName[imh], nuMass)
        histTmp = ROOT.TH2D(histname, histTitle, nbin_time, binning_time, nbins_Evis, Evismin, Evismax)
        h2d_etSpec.append( histTmp )
        h2d_etSpec[imh].GetXaxis().SetTitle("time [s]")
        h2d_etSpec[imh].GetYaxis().SetTitle("E_{vis} [MeV]")
        h2d_etSpec[imh].GetZaxis().SetTitle('%s event'%chaName[chaname])

        histname = 'hETBKG_mod%d_cha%d_mh%d'%(modelNum, chaname, imh)
        histTitle = '(T_{obs}, E_{vis}) for '
        histTitle += 'model%d, %s, %s, m_{#nu}%2.1f eV, with BKG'%(modelNum, chaName[chaname], MHName[imh], nuMass)
        histTmp = ROOT.TH2D(histname, histTitle, nbin_time, binning_time, nbins_Evis, Evismin, Evismax)
        h2d_etSpecBKG.append( histTmp )
        h2d_etSpecBKG[imh].GetXaxis().SetTitle("time [s]")
        h2d_etSpecBKG[imh].GetYaxis().SetTitle("E_{vis} [MeV]")
        h2d_etSpecBKG[imh].GetZaxis().SetTitle('%s event'%chaName[chaname])

        histname = 'hENuT_mod%d_cha%d_mh%d'%(modelNum, chaname, imh)
        histTitle = '(T_{obs}, E_{#nu}) for '
        histTitle += 'model%d, %s, %s, m_{#nu}%2.1f eV'%(modelNum, chaName[chaname], MHName[imh], nuMass)
        histTmp = ROOT.TH2D(histname, histTitle, nbin_time, binning_time, nbins_Ev, Evmin, Evmax)
        h2d_eNuTSpec.append( histTmp )
        h2d_eNuTSpec[imh].GetXaxis().SetTitle("time [s]")
        h2d_eNuTSpec[imh].GetYaxis().SetTitle("E_{v} [MeV]")
        h2d_eNuTSpec[imh].GetZaxis().SetTitle('%s event'%chaName[chaname])

        print('len(evisBins)-1', len(evisBins)-1)
        tmpSliceh2d = []
        for k in range(len(evisBins)-1):
            histname = 'hENuTSlice%d_mod%d_cha%d_mh%d'%(k, modelNum, chaname, imh)
            print(k, histname)
            histTitle = '(T_{obs}, E_{#nu}) for '
            histTitle += 'model%d, %s, %s, m_{#nu}%2.1f eV'%(modelNum, chaName[chaname], MHName[imh], nuMass)
            histTitle += ', E_{vis} in (%2.1f, %2.1f) MeV'%(evisBins[k], evisBins[k+1])
            histTmp = ROOT.TH2D(histname, histTitle, nbin_time, binning_time, nbins_Ev, Evmin, Evmax)
            #h2d_eNuTSpecSlice.append( histTmp )
            #h2d_eNuTSpecSlice[k].GetXaxis().SetTitle("time [s]")
            #h2d_eNuTSpecSlice[k].GetYaxis().SetTitle("E_{v} [MeV]")
            #h2d_eNuTSpecSlice[k].GetZaxis().SetTitle('%s event'%chaName[chaname])
            tmpSliceh2d.append(histTmp)
            tmpSliceh2d[k].GetXaxis().SetTitle("time [s]")
            tmpSliceh2d[k].GetYaxis().SetTitle("E_{v} [MeV]")
            tmpSliceh2d[k].GetZaxis().SetTitle('%s event'%chaName[chaname])
        h2d_eNuTSpecSlice.append(tmpSliceh2d)    

        # histname = 'hENuTCount_mod%d_cha%d_mh%d'%(modelNum, chaname, imh)
        # h2d_eNuTCount.append( ROOT.TH2D(histname, '', nbin_time, binning_time, nbins_Ev, Evmin, Evmax) )
        # histname = 'hETCount_mod%d_cha%d_mh%d'%(modelNum, chaname, imh)
        # h2d_etCount.append( ROOT.TH2D(histname, '', nbin_time, binning_time, nbins_Evis, Evismin, Evismax) )

    print('Nbin_Time=%d, Nbin_Evis=%d'%(nbin_time, nbins_Evis))

    # val = snDet.getEvisSpectrumAtTime(0.3, 0.3, nuType, MH)
    # print(val)

    # --- Interested ranges
    if chaname == 1:
        tBurstEnd = 1.0
    if chaname == 2:
        tBurstEnd = 1.0
    binBurstEnd = h2d_etSpec[0].GetXaxis().FindBin(tBurstEnd)
    if binBurstEnd>h2d_etSpec[0].GetNbinsX():
        binBurstEnd = h2d_etSpec[0].GetNbinsX()
    print('tBurstEnd, binBurstEnd: ', h2d_etSpec[0].GetXaxis().GetBinCenter(binBurstEnd), binBurstEnd)


    # ------------- Modify the PDFs due to different neutrino masses ---------------#
    Npar = 0
    # number of protons
    if chaname == 0 or chaname == 2:
        Npar = det.getNumberOfProton()
    # number of electrons
    if chaname == 1:
        Npar = det.getNumberOfProton()+6*det.getNumberOfCarbon()
    # number of carbon-12, 98.9% of carbons
    #if chaname == snDet.NCC or chaname == snDet.BCC or chaname == snDet.CNC:
    #    Npar = 0.989*peffLS.getNumberOfCarbon()
    print('Npar', Npar)

    timeTmp = array('d', [-999.]) # unit: s
    #me = 0.511; #MeV

    nCalc = 0
    Nabnormal = 0
    for ipt in range (binBurstEnd):
        itwidth = binning_time[ipt+1] - binning_time[ipt]
                
        #for k in range(nbins_Ev):
        #    EvTmp = Evmin + step_Ev * (k+0.5)
        for k in range(nbins_Evis):
            EvisTmp = Evismin + step_Evis * (k+0.5)
                            
            if nCalc%100000==0:
                print('nCalc: ', nCalc)
            nCalc += 1

            #EvTmpMin = 0.5*(EvisTmp+TMath::Sqrt(EvisTmp*EvisTmp+2*EvisTmp*me))
            #if chaname == 2:
            #    EvTmp = EvisTmp-0.51 + 1.3
            #    timeTmp[0] = binning_time[ipt]+0.5*itwidth 
            #    # timeTmp[0] is original neutrino time, 
            #    # after calling snFluenceDetAtTime, it will change due to non-zero neutrino mass
            #    # DeltaT = 5.14 * (mNu*mNu) * (100.0/Ev/Ev) * （snDist/10) ms
            #    fluence = modelSrc.snFluenceDetAtTime(timeTmp, nuMass, EvTmp, nuType, MH)
            #    xs = peffLS.totalXS(EvTmp, nuType)
            #    if fluence > 0. and xs >0.:
            #        tBin = dummyT.FindBin(timeTmp[0])
            #        tFactor = itwidth/dummyT.GetBinWidth(tBin)

            #        evisFactor = 1
            #        h2d_etSpec[MH].Fill(timeTmp[0], EvisTmp, Npar*fluence*xs*evisFactor*tFactor)

            #        evisFactor = 1
            #        h2d_eNuTSpec[MH].Fill(timeTmp[0], EvTmp, Npar*fluence*xs*evisFactor*tFactor)
            #        # h2d_eNuTCount[MH].Fill(timeTmp[0], EvTmp)
            #        # h2d_etCount[MH].Fill(timeTmp[0], EvisTmp)
            #        # h2d_etSpec[MH].Fill(timeTmp[0], EvisTmp, Npar*fluence*xs*itwidth*step_Evis)
            #        # h2d_eNuTSpec[MH].Fill(timeTmp[0], EvTmp, Npar*fluence*xs*itwidth*step_Evis)


            #if chaname == 1:   # Nue
            #    t0 = binning_time[ipt]+0.5*itwidth
            #    fluence, tCorr = SNGar.oneSNFluenceDetTimeShift(t0, nuMass, EvTmp, nuType, MH)
            #    xs = nuexs.totalXS(EvTmp, nuType)
            #    print(ipt, t0, k, EvTmp, xs, fluence)
            #    if fluence > 0 and xs > 0:
            #        tBin = dummyT.FindBin(timeTmp[0])
            #        tFactor = itwidth/dummyT.GetBinWidth(tBin)

            #        evisFactor = 1
            #        h2d_eNuTSpec[MH].Fill(tCorr, EvTmp, Npar*fluence*xs*evisFactor*tFactor)

    


            if chaname == 1:
                for j in range(nbins_Ev):
                    EvTmp = Evmin + (j+0.5) * step_Ev
                    #timeTmp[0] = binning_time[ipt]+0.5*itwidth
                    t0 = binning_time[ipt]+0.5*itwidth
                    # timeTmp[0] is original neutrino time, 
                    # after calling snFluenceDetAtTime, it will change due to non-zero neutrino mass
                    #fluence = modelSrc.snFluenceDetAtTime(timeTmp, nuMass, EvTmp, nuType, MH)
                    fluence, tCorr = SNGar.oneSNFluenceDetTimeShift(t0, nuMass, EvTmp, nuType, MH)
                    #dxs = peffLS.differentialXS(EvTmp, EvisTmp, nuType)
                    dxs = nuexs.differentialXS(EvTmp, EvisTmp, nuType)
                    #print(nCalc, t0, tCorr, EvTmp, EvisTmp, dxs, fluence)
                    if fluence > 0. and dxs >0.:
                        tBin = dummyT.FindBin(timeTmp[0])
                        tFactor = itwidth/dummyT.GetBinWidth(tBin)

                        evisFactor = step_Ev/step_Evis
                        h2d_etSpec[MH].Fill(tCorr, EvisTmp, Npar*fluence*dxs*evisFactor*tFactor)

                        evisFactor = step_Evis/step_Ev
                        h2d_eNuTSpec[MH].Fill(tCorr, EvTmp, Npar*fluence*dxs*evisFactor*tFactor)

                    #    if EvisTmp<evisBins[len(evisBins)-1] and EvisTmp>evisBins[0]:
                    #        eBin = dummyEvis.FindBin(EvisTmp)
                    #        h2d_eNuTSpecSlice[MH][eBin-1].Fill(timeTmp[0], EvTmp, Npar*fluence*dxs*evisFactor*tFactor)

                            # if h2d_eNuTSpecSlice[len(evisBins)-2].GetEntries()>Nabnormal:
                            #     print('test1: ', eBin, h2d_eNuTSpecSlice[10].GetEntries(), timeTmp[0], EvisTmp, EvTmp)
                            #     Nabnormal = Nabnormal + 1

                        #### Save different eNuTSpec according to Evis
                        #h2d_eNuTCount[MH].Fill(timeTmp[0], EvTmp, step_Evis)
                        #h2d_etCount[MH].Fill(timeTmp[0], EvisTmp)
                # if h2d_eNuTSpecSlice[10].GetEntries()>0:
                #     print('error')
    # -------------------------------------------------------------------------


    fout.cd()

    tinfo.Write()
    tmpBKG = 60/86400./15
    for imh in range(nMH):
        if imh!=MH:
            continue
        # h2d_eNuTCount[imh].Write()
        # h2d_etCount[imh].Write()
        # normalize
        # if chaname==snDet.IBD: 
        #     for i in range(1, h2d_eNuTCount[imh].GetNcells()+1):
        #         if h2d_eNuTCount[imh].GetBinContent(i)!=0:
        #             val = h2d_eNuTSpec[imh].GetBinContent(i)
        #             count = h2d_eNuTCount[imh].GetBinContent(i)
        #             h2d_eNuTSpec[imh].SetBinContent(i, val/count)
        h2d_eNuTSpec[imh].Write()

        # if chaname==snDet.IBD: 
        #     for i in range(1, h2d_etCount[imh].GetNcells()+1):
        #         if h2d_etCount[imh].GetBinContent(i)!=0:
        #             val = h2d_etSpec[imh].GetBinContent(i)
        #             count = h2d_etCount[imh].GetBinContent(i)
        #             h2d_etSpec[imh].SetBinContent(i, val/count)
        h2d_etSpec[imh].Write()

        #for k in range(len(evisBins)-1):
        #    h2d_eNuTSpecSlice[MH][k].Write()

        # Set background from reactors and others
        Xbins = h2d_etSpecBKG[imh].GetNbinsX()
        Ybins = h2d_etSpecBKG[imh].GetNbinsY()
        for i in range(1, Xbins+1):
            Xwidth = h2d_etSpecBKG[imh].GetXaxis().GetBinWidth(i)

            for j in range(1, Ybins+1):
                if h2d_etSpecBKG[imh].GetYaxis().GetBinCenter(j)<1.0:
                    continue

                Ywidth = h2d_etSpecBKG[imh].GetYaxis().GetBinWidth(j)
                #val = h2d_etSpec[imh].GetBinContent(i, j)

                h2d_etSpecBKG[imh].SetBinContent(i, j, tmpBKG)
        
        h2d_etSpecBKG[imh].Write()

        print('Integral(): imh, et, eNuT', imh, h2d_etSpec[imh].Integral("width"), h2d_eNuTSpec[imh].Integral("width"))

    fout.Close()

    #fig2 = plt.figure(figsize=(5,5))
    #ax = fig2.add_subplot()
    ##ec = plt.hist(relTobsNuebar, 50, histtype='stepfilled', facecolor='g', range=(-1,9))
    #ec = plt.hist(nuNum, 50, histtype='step', facecolor='b', range=(0,5))
    #ax.set_xlabel('E_{vis} (MeV)')
    #plt.show()
    #print("test")

