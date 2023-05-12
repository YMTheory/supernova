import numpy as np
import argparse
import ROOT
import os
import pandas as pd
import matplotlib as mpl

#import multiprocessing
#from multiprocessing import cpu_count

from channel_analyser import channel
import scanner

import sys
sys.path.append("/junofs/users/miaoyu/supernova/simulation/toyMC/script/")
import stat_tool

import matplotlib.pyplot as plt
import plotConfig


def readCSV(filename, start, end):
    df1 = pd.read_csv(filename)
    Tbest = df1["Tbest"][start:end]
    locMin = df1["locMin"][start:end]
    return Tbest, locMin


def loadFit(MO, mass, evtID):
    csvfilepath = os.getenv("CSVFILEPATH")
    filename = f"{csvfilepath}Garching82703_10.0kpc_20.0kton_{MO}_pESeESIBD_nuMass{mass:.1f}eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_noC14_tot_unbinnedData_unbinnedNLL_PosiToyDataTobs2D_JUNO.csv"
    tmp_Tbest, tmp_locMin = readCSV(filename, evtID, evtID+1)
    return tmp_Tbest, tmp_locMin


if __name__ == "__main__" :
    """
    Feasible diagnosis tool, mainly to understand the absolute mass analysis.
    """

    MO          = "NO"
    model       = "Garching"
    modelNo     = 82703
    Ethr        = 0.10
    fitTmin     = -0.02
    fitTmax     = 0.02
    nuMass      = 0.0
    fileNo      = 0
    dist        = 10.
    target      = 1.
    exp         = "JUNO"
    startevt    = 0
    endevt      = 100
    doFit       = True
    fitDim      = 2
    asimov      = False
    eES         = True
    IBD         = True
    pES         = True
    C14         = False
    C14level    ="low"
    plot        = True

    parser = argparse.ArgumentParser(description='Arguments of SNNu analyser.')
    parser.add_argument('--model',      type=str,       default="Garching",     help="Model name (option: Garching, Burrows).")
    parser.add_argument('--modelNo',    type=int,       default=82703,          help="Model number.")
    parser.add_argument('--MO',         type=str,       default="NO",           help="Mass ordering for the dataset.")
    parser.add_argument("--start",      type=int,       default=0,              help="Start event for fit.")
    parser.add_argument("--end",        type=int,       default=0,              help="End event for fit.")
    parser.add_argument("--fileNo",     type=int,       default=0,              help="File No.")
    parser.add_argument("--doFit",      dest="fit",     action="store_true",    help="enable fitting.")
    parser.add_argument("--no-doFit",   dest="fit",     action="store_false",   help="disable fitting.")
    parser.add_argument("--doPlot",      dest="plot",     action="store_true",    help="enable plotting.")
    parser.add_argument("--no-doPlot",   dest="plot",     action="store_false",   help="disable plotting.")
    parser.add_argument("--Asimov",     dest="asimov",  action="store_true",    help="enable asimov dataset.")
    parser.add_argument("--no-Asimov",  dest="asimov",  action="store_false",   help="disable asimov dataset.")
    parser.add_argument("--C14",        dest="C14",     action="store_true",    help="enable C14 background.")
    parser.add_argument("--no-C14",     dest="C14",     action="store_false",   help="disable C14 background.")
    parser.add_argument("--fitDim",     type=int,       default=2,              help="Fitting dimensions (1 for time only, 2 to time combining energy.)")
    parser.add_argument("--C14level",   type=str,       default="low",          help="C14 concentration level (low or high).")
    parser.add_argument("--Ethr",       type=float,     default=0.1,            help="Detection threshold MeV. ")
    parser.add_argument("--nuMass",     type=float,     default=0.0,            help="Neutrino mass in fit PDFs.")
    parser.add_argument("--fitTmin",    type=float,     default=-0.02,          help="Minimum fit time [s].")
    parser.add_argument("--fitTmax",    type=float,     default=0.02,           help="Maximum fit time [s].")
    parser.add_argument('--eES',        dest='eES',     action="store_true",    help="enable eES")
    parser.add_argument('--no-eES',     dest='eES',     action="store_false",   help="disable eES")
    parser.add_argument('--IBD',        dest='IBD',     action="store_true",    help="enable IBD")
    parser.add_argument('--no-IBD',     dest='IBD',     action="store_false",   help="disable IBD")
    parser.add_argument('--pES',        dest='pES',     action="store_true",    help="enable pES")
    parser.add_argument('--no-pES',     dest='pES',     action="store_false",   help="disable pES")
    parser.add_argument('--dist',       type=float,     default=10.,            help="distance of CCSNe.")
    parser.add_argument('--target',     type=float,     default=1.,             help="normalized target mass, JUNO 20kton scaled as 1.")
    args = parser.parse_args()

    model       = args.model
    modelNo     = args.modelNo
    MO          = args.MO
    startevt    = args.start
    endevt      = args.end
    doFit       = args.fit
    fileNo      = args.fileNo
    fitDim      = args.fitDim
    nuMass      = args.nuMass
    asimov      = args.asimov
    eES         = args.eES
    IBD         = args.IBD
    pES         = args.pES
    Ethr        = args.Ethr
    fitTmin     = args.fitTmin
    fitTmax     = args.fitTmax
    C14         = args.C14
    C14level    = args.C14level
    dist        = args.dist
    target      = args.target
    plot        = args.plot

    scale = target * 10.0 * 10.0 / (dist * dist)  ## due to the PDFs are generated based on 10 kpc CCSNe and 20kton LS...
    bkgscale = target

    if Ethr > 0.2:
        Eibd = Ethr
    else:
        Eibd = 0.2


    channels = {}
    if pES:
        channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fitEmax=5, fileNo=fileNo, dist=dist, exp=exp)
        channels["pES"].setC14rate(0)
    if IBD:
        channels["IBD"] = channel("IBD", MO, model, modelNo, Eibd, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
        channels["IBD"].setC14rate(0)
    if eES:
        channels["eES"] = channel("eES", MO, model, modelNo, Eibd, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
        channels["eES"].setC14rate(0)

    for ich, cha in enumerate(channels.values()):

        cha.setNevtPerFile(1e5)
        cha.setStartEvtId(startevt)
        cha.setEndEvtId(endevt)

        cha.setScale(scale)
        cha.setBkgScale(bkgscale)

        # 1D PDF setting:

        if cha.name == "eES" or cha.name == "IBD" or cha.name == "pES":
            pdffile1dpath = os.getenv("PDF1DFILEPATH")
            cha.setDataPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_{cha.MH}_10kpc_{cha.name}_nuMass0.0eV_TEobs2dPDFintegral_JUNO.root")
            cha._load_datapdf1D()
            cha.setNOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDFintegral_JUNO.root")
            cha.setIOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDFintegral_JUNO.root")
            cha._load_pdf1D()

            if fitDim == 2:
                pdffile2dpath = os.getenv("PDF2DFILEPATH")
                print(f"2D PDF file path: {pdffile2dpath}")
                cha.setDataPdfFile2D(f"{pdffile2dpath}{model}{modelNo}_nuePDF_{cha.MH}_10kpc_{cha.name}_nuMass0.0eV_TEobs2dPDF_JUNO.root")
                cha._load_datapdf2D()
                cha.setNOPdfFile2D(f"{pdffile2dpath}{model}{modelNo}_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_JUNO.root")
                cha.setIOPdfFile2D(f"{pdffile2dpath}{model}{modelNo}_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_JUNO.root")
                cha._load_pdf2D()

        # Set Data Files
        if not asimov:
            if fitDim == 1:
                cha.setDataFile1D(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data1d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_T1D.root")
                cha._load_data1D()
            if fitDim == 2:
                if cha.name == "pES":
                    if not C14:
                        cha.setDataFile2D(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_noC14.root")
                    else:
                        if C14level == "low":
                            cha.setDataFile2D(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_C14low.root")
                        elif C14level == "high":
                            cha.setDataFile2D(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_C14high.root")
                else:
                    datapath = os.getenv("DATA2DFILEPATH")
                    cha.setDataFile2D(f"{datapath}Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_JUNO.root")

                cha._load_data2D()    # could get 1D or 2D dataset from the fitting requirement.

        FITTING_EVENT_NUM =  cha.getNevtPerFile() # the sample number to run...
    

    
    if plot:
        FIGSIZE = 4
        colors0 = mpl.cm.viridis(np.linspace(0.1,1, 21))
        fig, ax = plt.subplots(1, len(channels), figsize=(len(channels)*FIGSIZE+FIGSIZE, FIGSIZE) )
        for ich, cha in enumerate(channels.values()):
            for ievt in range(startevt, endevt, 1):
                # load and draw data
                oneT, oneE = cha.get_one_event2D(ievt)
                conts, edges = np.histogram(np.array(oneT), bins=40, range=(-0.02, 0.02))
                dx, dy, dyerr = [], [], []
                for idd in range(len(conts)):
                    if conts[idd] > 0:
                        dx.append((edges[idd]+edges[idd+1])/2.)
                        dy.append(conts[idd])
                        dyerr.append(np.sqrt(conts[idd]))
                ax[ich].errorbar(dx, dy, yerr=dyerr, fmt="o", ms=5, fillstyle="none", color="black")
                # load and draw PDFs
                for im, mass in enumerate(np.arange(0, 2.0, 0.1)):
                    dTbest, _ = loadFit("NO", mass, startevt)
                    cha.setNOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_NO_10kpc_{cha.name}_nuMass{mass:.1f}eV_TEobs2dPDFintegral_JUNO.root")
                    cha.setIOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_IO_10kpc_{cha.name}_nuMass{mass:.1f}eV_TEobs2dPDFintegral_JUNO.root")
                    x = np.arange(-0.02, 0.02, 0.001)
                    y = np.zeros(len(x))
                    for ix, ex in enumerate(x):
                        y[ix] = cha._pdfNO_func(ex+dTbest)
                    ax[ich].plot(x, y/1000., ":", lw=1.5, label=f"{mass:.1f}eV", color=colors0[im])
            plotConfig.setaxis(ax[ich], xlabel="Post-bounce time [ms]", ylabel="Counts per ms", title=cha.name, lg=True, lgsize=10, ncol=2)
        #plt.tight_layout()
        plt.show()

                    

    
