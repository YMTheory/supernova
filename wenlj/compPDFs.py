#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: Liangjian Wen ---> wenlj@ihep.ac.cn
# Created Time : Sun Nov 1
# Compare different (time, Ev) PDFs
"""

import numpy as np
import matplotlib.pyplot as plt
import ROOT
#ROOT.gSystem.Load("/mnt/c/Users/LiangjianWen/Documents/JUNO/Physics/SNsim/simulation/lib/libSNsim.so")
from ROOT import TH1D, TH2D, TCanvas
import scipy.stats as st
from scipy.stats import rv_continuous
from array import array
import sys

def findBin(hist2d, x1, x2, y1, y2):
    binx1 = hist2d.GetXaxis().FindBin(x1)
    binx2 = hist2d.GetXaxis().FindBin(x2)
    biny1 = hist2d.GetYaxis().FindBin(y1)
    biny2 = hist2d.GetYaxis().FindBin(y2)
    return (binx1, binx2, biny1, biny2)

if __name__ == "__main__":

    # file = ROOT.TFile()

    # fig2 = plt.figure(figsize=(5,5))
    # ax = fig2.add_subplot()
    # ec = plt.hist(relTobsNuebar, 50, histtype='stepfilled', facecolor='g', range=(-1,9))
    # ec = plt.hist(nuNum, 50, histtype='step', facecolor='b', range=(0,5))
    # ax.set_xlabel('E_{vis} (MeV)')
    # plt.show()
    # print("test")
    ROOT.gROOT.SetBatch(True)

    if len(sys.argv) < 7:
        print("Usage: %s [modelNum] [channel] [nuType] [nuMass1] [nuMass2] [dist]"  % (sys.argv[0]))
        print("modelNum: check the SNsim/simulation/data directory")
        print("channel : enum {NuP, NuE, IBD, NCC, BCC, CNC}")
        print("nuType  : -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x")
        print("nuMass1 : neutrino mass that used to generate dataset, in eV unit")
        print("nuMass2 : neutrino mass that used in PDF to fit the dataset, in eV unit")
        print("dist    : SN distance in kpc unit")
        print("example: python compPDFs.py 82503 1 0 0.0 0.2 0.4")
        sys.exit(1)

    ##################################
    # Input variables
    ROOT.RooFit.PrintLevel(1)
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
    
    nuMass1 = float( sys.argv[4] )# unit: eV
    print('nuMass for producing dataset : ', nuMass1, ' eV')

    nuMass2 = float( sys.argv[5] )# unit: eV
    print('nuMass for producing PDF: ', nuMass2, ' eV')

    dist = float( sys.argv[6] )# unit: kpc
    print('dist: ', dist, ' kpc')

    MH = 1
    ##################################
    can = ROOT.TCanvas("can", "can", 1200, 800)
    can.Divide(1,2)
    can.cd(1).SetTopMargin(0.04)
    can.cd(2).SetTopMargin(0.04)
    can.cd(1).SetLeftMargin(0.12)
    can.cd(2).SetLeftMargin(0.12)
    #pdfilename = 'compEtruePDFs_data%2.1feV_pdf%2.1feV_%2.1fkpc.pdf'%(nuMass1,nuMass2,dist)
    pdfilename = './spectra/compEPDFs_data%2.1feV_pdf%2.1feV_%2.1fkpc_MH%d_2.pdf'%\
                 (nuMass1,nuMass2,dist,MH)
    can.Print(pdfilename+'[')
    # Obtain histograms from the root files
    prefix = "./etSpec/fineSpec/TEvisPDF_"
    prefix += 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)
    filename = prefix+'_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25.root'%(nuMass1, dist)
    print(filename)
    dataFile1  = ROOT.TFile(filename, 'read')
    #inputHist1 = dataFile1.Get("hENuT_mod%d_cha%d_mh0"%(modelNum, chaname))
    inputHist1 = dataFile1.Get("hET_mod%d_cha%d_mh%d"%  (modelNum, chaname, MH))
    hENuT1     = dataFile1.Get("hENuT_mod%d_cha%d_mh%d"%(modelNum, chaname, MH))

    filename = prefix+'_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25.root'%(nuMass2, dist)
    print(filename)
    dataFile2  = ROOT.TFile(filename, 'read')
    #inputHist2 = dataFile2.Get("hENuT_mod%d_cha%d_mh0"%(modelNum, chaname))
    inputHist2 = dataFile2.Get("hET_mod%d_cha%d_mh%d"%  (modelNum, chaname, MH))
    hENuT2     = dataFile2.Get("hENuT_mod%d_cha%d_mh%d"%(modelNum, chaname, MH))

    nuMass3 = 0.6
    filename = prefix+'_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25.root'%(nuMass3, dist)
    print(filename)
    dataFile3  = ROOT.TFile(filename, 'read')
    #inputHist2 = dataFile2.Get("hENuT_mod%d_cha%d_mh0"%(modelNum, chaname))
    inputHist3 = dataFile3.Get("hET_mod%d_cha%d_mh%d"%  (modelNum, chaname, MH))
    hENuT3     = dataFile3.Get("hENuT_mod%d_cha%d_mh%d"%(modelNum, chaname, MH))

    ###########################################################
#    ##### Plot the 2D histogram
#    #fig = plt.figure(figsize=(5,5))
#    #dataENuT = np.asarray(hENuT2)
#    #print(dataENuT[1])
#
#    #ax = fig.add_subplot()
#    #ec = plt.hist2d(dataENuT[:,0], dataENuT[:,1], cmap='jet')
#    #ax.set_ylabel('$E_{\nu}$ (MeV)')
#    #plt.show()
#    #print("test")
#    #plt.savefig('./spectra/fineSpec/plotPDF_%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d.pdf'%(nuMass2, dist, Ethr, group), dpi=400, bbox_inches='tight')

    can2 = ROOT.TCanvas("can2", "can2", 1200, 600)
    can2.SetTopMargin(0.04)
    can2.SetRightMargin(0.12)
    can2.SetLeftMargin(0.08)
    can2.SetBottomMargin(0.13)
    ROOT.gPad.SetLogz(1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(55)

    ############################
    if MH==0:
        hENuT1.GetXaxis().SetRangeUser(-0.025, 0.025)
    if MH==1:
        hENuT1.GetXaxis().SetRangeUser(-0.025, 0.05)
    if MH==2:
        hENuT1.GetXaxis().SetRangeUser(-0.025, 0.05)
    hENuT1.GetYaxis().SetRangeUser(0.0, 40.0)
    hENuT1.GetYaxis().SetTitleOffset(0.6)
    hENuT1.GetYaxis().SetNdivisions(505)
    hENuT1.Draw('colz')
    pt = ROOT.TPaveText(0.15, 0.8, 0.35, 0.95, 'NB NDC')
    pt.SetFillColor(0)
    pt.AddText('m_{#nu} = %2.1f eV'%nuMass1)
    pt.SetTextSize(0.06)
    pt.Draw()
    picName = './spectra/fineSpec/plotPDF_cha%d%s_%2.1feV_%2.1fkpc_MH%d'\
              %(chaname, chaName[chaname], nuMass1, dist, MH)
    can2.Print(picName+'.pdf')
    can2.Print(picName+'.png')

    ############################
    if MH==0:
        hENuT2.GetXaxis().SetRangeUser(-0.025, 0.025)
    if MH==1:
        hENuT2.GetXaxis().SetRangeUser(-0.025, 0.05)
    if MH==2:
        hENuT2.GetXaxis().SetRangeUser(-0.025, 0.05)
    hENuT2.GetYaxis().SetRangeUser(0.0, 40.0)
    hENuT2.GetYaxis().SetTitleOffset(0.6)
    hENuT2.GetYaxis().SetNdivisions(505)
    hENuT2.Draw('colz')

    pt = ROOT.TPaveText(0.15, 0.8, 0.35, 0.95, 'NB NDC')
    pt.SetFillColor(0)
    pt.AddText('m_{#nu} = %2.1f eV'%nuMass2)
    pt.SetTextSize(0.06)
    pt.Draw()
    picName = './spectra/fineSpec/plotPDF_cha%d%s_%2.1feV_%2.1fkpc_MH%d'\
              %(chaname, chaName[chaname], nuMass2, dist, MH)
    can2.Print(picName+'.pdf')
    can2.Print(picName+'.png')

    ############################
    if MH==0:
        hENuT3.GetXaxis().SetRangeUser(-0.025, 0.025)
    if MH==1:
        hENuT3.GetXaxis().SetRangeUser(-0.025, 0.05)
    if MH==2:
        hENuT3.GetXaxis().SetRangeUser(-0.025, 0.05)
    hENuT3.GetYaxis().SetRangeUser(0.0, 40.0)
    hENuT3.GetYaxis().SetTitleOffset(0.6)
    hENuT3.GetYaxis().SetNdivisions(505)
    hENuT3.Draw('colz')

    pt = ROOT.TPaveText(0.15, 0.8, 0.35, 0.95, 'NB NDC')
    pt.SetFillColor(0)
    pt.AddText('m_{#nu} = %2.1f eV'%nuMass3)
    pt.SetTextSize(0.06)
    pt.Draw()
    picName = './spectra/fineSpec/plotPDF_cha%d%s_%2.1feV_%2.1fkpc_MH%d'\
              %(chaname, chaName[chaname], nuMass3, dist, MH)
    can2.Print(picName+'.pdf')
    can2.Print(picName+'.png')

    ############################
    can3 = ROOT.TCanvas("can3", "can3", 1200, 600)
    can3.SetTopMargin(0.04)
    can3.SetRightMargin(0.03)
    can3.SetLeftMargin(0.1)
    can3.SetBottomMargin(0.13)

    Tmin, Tmax = -0.05, 0.05
    if MH==1 or MH==2:
        Tmax = 0.1
    Emin, Emax = 0.1, 10.0
    if chaname==2:
        Emin, Emax = 1.02, 10+1.02
        #Emin, Emax = 1.02, 25+1.02

    binx1, binx2, biny1, biny2 = findBin(inputHist1, Tmin, Tmax, Emin, Emax)
    projX1 = inputHist1.ProjectionX('_px1', biny1, biny2, 'E')
    projX1.SetLineColor(4)
    histTitle = 'PDF comparison, blue %2.1feV, red %2.1feV,%2.1fkpc'%(nuMass1, nuMass2, dist)
    histTitle += ', E_{vis} in (%2.1f MeV, %2.1f MeV)'%(Emin, Emax)
    #histTitle += ', E_{#nu} in (%2.1f MeV, %2.1f MeV)'%(Emin, Emax)
    projX1.SetTitle(histTitle)
    projX1.GetXaxis().SetRangeUser(Tmin, Tmax)
    projX2 = inputHist2.ProjectionX('_px2', biny1, biny2, 'E')
    projX2.GetXaxis().SetRangeUser(Tmin, Tmax)
    projX2.SetLineColor(2)
    projX2.SetLineStyle(2)

    projX3 = inputHist3.ProjectionX('_px3', biny1, biny2, 'E')
    projX3.GetXaxis().SetRangeUser(Tmin, Tmax)
    projX3.SetLineColor(1)
    projX3.SetLineStyle(7)

    projX1.GetYaxis().SetTitle('Events [ms^{-1}]')
    projX1.GetYaxis().SetTitleOffset(0.8)

    scale = 5e-2 * 5e-2 # 0.05 ms per time bin, 0.05 MeV per Energy bin
    projX1.Scale(scale)
    projX2.Scale(scale)
    projX3.Scale(scale)

    projX1.Draw("hist")
    projX2.Draw("same hist")
    projX3.Draw("same hist")
    projX1.GetXaxis().SetTitleOffset(1.0)

    leg = ROOT.TLegend(0.15, 0.6, 0.4, 0.85)
    pt = leg.AddEntry(projX1, 'm_{#nu} = %2.1f eV'%nuMass1, 'L')
    pt.SetTextColor(4)
    pt.SetTextFont(132)
    pt = leg.AddEntry(projX3, 'm_{#nu} = %2.1f eV'%nuMass3, 'L')
    pt.SetTextColor(1)
    pt.SetTextFont(132)
    pt = leg.AddEntry(projX2, 'm_{#nu} = %2.1f eV'%nuMass2, 'L')
    pt.SetTextColor(2)
    pt.SetTextFont(132)
    leg.Draw()
    picName = './spectra/fineSpec/compPDF_cha%d%s_%2.1feV_m%2.1feV_%2.1fkpc_MH%d_Emax%3.1f'\
              %(chaname, chaName[chaname], nuMass1, nuMass2, dist, MH, Emax)
    can3.Print(picName+'.pdf')
    can3.Print(picName+'.png')

    ###########################################################

#    Tmin, Tmax = -0.05, 0.05
#
#    # Evis
#    #Ebins = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 
#    #          4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
#    Ebins = [0.1, 1.0]
#
#    # Ebins = [0., 1.0, 2.0, 3.0, 
#    #           4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10, 11, 12, 13]
#
#
#    for k in range(len(Ebins)-1):
#        Emin, Emax = Ebins[k], Ebins[k+1]
#
#        binx1, binx2, biny1, biny2 = findBin(inputHist1, Tmin, Tmax, Emin, Emax)
#        # projY1 = inputHist1.ProjectionY('_py1', binx1, binx2, '')
#        # projY1.SetLineColor(1)
#        # can.cd(1)
#        # projY1.SetTitle('PDF comparison, black %2.1feV, red %2.1feV'%(nuMass1, nuMass2))
#        # projY1.Draw()
#        # projY2 = inputHist2.ProjectionY('_py2', binx1, binx2, '')
#        # projY2.SetLineColor(2)
#        # projY2.Draw('same')
#
#        can.cd(1)
#        projX1 = inputHist1.ProjectionX('_px1_%d'%(k), biny1, biny2, 'E')
#        projX1.SetLineColor(4)
#        histTitle = 'PDF comparison, blue %2.1feV, red %2.1feV,%2.1fkpc'%(nuMass1, nuMass2, dist)
#        histTitle += ', E_{vis} in (%2.1f MeV, %2.1f MeV)'%(Emin, Emax)
#        #histTitle += ', E_{#nu} in (%2.1f MeV, %2.1f MeV)'%(Emin, Emax)
#        projX1.SetTitle(histTitle)
#        projX1.SetTitleOffset(0.01)
#        projX1.GetXaxis().SetRangeUser(Tmin, Tmax)
#        projX2 = inputHist2.ProjectionX('_px2_%d'%(k), biny1, biny2, 'E')
#        projX2.GetXaxis().SetRangeUser(Tmin, Tmax)
#        projX2.SetLineColor(2)
#        projX2.SetTitle(histTitle)
#        projX1.Draw("hist")
#        projX2.Draw("same hist")
#
#        can.cd(2)
#        hDiff = projX1.Clone('hDiff_%d'%(k))
#        hDiff.Reset('ICE')
#        for j in range(1, hDiff.GetNbinsX()+1):
#            count1 = projX1.GetBinContent(j)
#            count2 = projX2.GetBinContent(j)
#            mean = (count2+count1)/2
#            if mean!=0:
#                ratio = (count2-count1)/mean
#                hDiff.SetBinContent(j, ratio)
#        hDiff.GetYaxis().SetRangeUser(-0.4, 0.4)
#        hDiff.GetYaxis().SetTitle('diff in (%2.1f, %2.1f) MeV'%(Emin, Emax))
#        hDiff.GetYaxis().SetTitleOffset(0.7)
#        hDiff.Draw("hist")
#        
#        can.Print(pdfilename)
#
#    # dataFile1.cd()
#    # projY1.GetXaxis().SetRangeUser(Emin, Emax)
#    # can.cd(1)
#    # projY1.Draw()
#    # dataFile2.cd()
#    # projY2.Draw('same')
#    
#    # dataFile1.cd()
#    # projX1.GetXaxis().SetRangeUser(Tmin, Tmax)
#    # can.cd(2)
#    # projX1.Draw()
#    # dataFile2.cd()
#    # projX2.Draw('same')
#
#    #can.Print(pdfilename)
#
#    # 
#    Tbins = [-0.01, 0.002, 0.01]
#    Emin, Emax = 0.1, 4.0
#    for k in range(len(Tbins)-1):
#        Tmin, Tmax = Tbins[k], Tbins[k+1]
#
#        binx1, binx2, biny1, biny2 = findBin(inputHist1, Tmin, Tmax, Emin, Emax)
#        # projY1 = inputHist1.ProjectionY('_py1', binx1, binx2, '')
#        # projY1.SetLineColor(1)
#        # can.cd(1)
#        # projY1.SetTitle('PDF comparison, black %2.1feV, red %2.1feV'%(nuMass1, nuMass2))
#        # projY1.Draw()
#        # projY2 = inputHist2.ProjectionY('_py2', binx1, binx2, '')
#        # projY2.SetLineColor(2)
#        # projY2.Draw('same')
#
#        can.cd(1)
#        projY1 = inputHist1.ProjectionY('_py1_%d'%(k), binx1, binx2, 'E')
#        projY1.SetLineColor(4)
#        histTitle = 'PDF comparison, blue %2.1feV, red %2.1feV,%2.1fkpc'%(nuMass1, nuMass2, dist)
#        histTitle += ', T in (%2.1f s, %2.1f s)'%(Tmin, Tmax)
#        #histTitle += ', E_{#nu} in (%2.1f MeV, %2.1f MeV)'%(Emin, Emax)
#        projY1.SetTitle(histTitle)
#        projY1.GetXaxis().SetRangeUser(Tmin, Tmax)
#        projY2 = inputHist2.ProjectionY('_py2_%d'%(k), binx1, binx2, 'E')
#        projY2.GetXaxis().SetRangeUser(Tmin, Tmax)
#        projY2.SetLineColor(2)
#        scale = 1./projY1.Integral()
#        projY1.Scale(scale)
#        scale = 1./projY2.Integral()
#        projY2.Scale(scale)
#        projY1.GetYaxis().SetRangeUser(0., 0.014)
#        projY1.Draw("hist")
#        projY2.Draw("same hist")
#
#        can.cd(2)
#        hDiff = projY1.Clone('hDiff_%d'%(k))
#        hDiff.Reset('ICE')
#        for j in range(1, hDiff.GetNbinsX()+1):
#            count1 = projY1.GetBinContent(j)
#            count2 = projY2.GetBinContent(j)
#            mean = (count2+count1)/2
#            if mean!=0:
#                ratio = (count2-count1)/mean
#                hDiff.SetBinContent(j, ratio)
#        hDiff.GetYaxis().SetRangeUser(-0.05, 0.05)
#        hDiff.GetYaxis().SetTitle('diff in (%2.1f, %2.1f) MeV'%(Tmin, Tmax))
#        hDiff.GetYaxis().SetTitleOffset(0.9)
#        hDiff.Draw("hist")
#        
#        can.Print(pdfilename)
#
#    can.Print(pdfilename+']')
#
