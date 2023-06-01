import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm
import ROOT
import warnings
warnings.filterwarnings("ignore")
import logging

import multiprocessing
from multiprocessing import cpu_count
import time

from channel_analyser import channel
import scanner

import os
import sys
sys.path.append("/junofs/users/miaoyu/supernova/simulation/toyMC/script/")
import stat_tool


def fit_absolute_mass():

    """
    This script is used to do combining fitting: neutrino mass ordering and absolute neutrion mass.
    - The argument --useMass is to specify if it is the absolute neutrino mass fitting or mass ordering fitting.
    For a specific type of fitting, there is two methods: asimov dataset test and toyMC fitting.
    - The argument --Asimov is to specify which method to use.

    Detector configuration:
    - Argument --target: scale factor of target mass 20 kton, e.g. target = 0.5 means target mass = 20 * 0.5 = 10 kton, C14 event is also scaled with this factor.
    - Argument --dist: distance of CCSNe.


    Background configuration:
    - Argument --C14 to specify whether considering C14 background;
    - Argument --C14level to specify the C14 concentration level: currently 

    Fitting configuration:
    - Argument --fitDim to specify the fitting dimensions, 1 for all time spectral fit, 2 for only pES 2D fit and eES/IBD still 1D fit.
    """

    logging.basicConfig(filename='example.log',level=logging.DEBUG)

    print("==================================\n")
    print(f"Total CPU number is {cpu_count()}")
    print("==================================\n")

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
    output      = True
    fitDim      = 2
    asimov      = False
    eES         = True
    IBD         = True
    pES         = True
    C14         = False
    C14level    ="low"
    scanning    = False

    parser = argparse.ArgumentParser(description='Arguments of SNNu analyser.')
    parser.add_argument('--model',      type=str,       default="Garching",     help="Model name (option: Garching, Burrows).")
    parser.add_argument('--modelNo',    type=int,       default=82703,          help="Model number.")
    parser.add_argument('--MO',         type=str,       default="NO",           help="Mass ordering for the dataset.")
    parser.add_argument("--start",      type=int,       default=0,              help="Start event for fit.")
    parser.add_argument("--end",        type=int,       default=0,              help="End event for fit.")
    parser.add_argument("--fileNo",     type=int,       default=0,              help="File No.")
    parser.add_argument("--doFit",      dest="fit",     action="store_true",    help="enable fitting.")
    parser.add_argument("--no-doFit",   dest="fit",     action="store_false",   help="disable fitting.")
    parser.add_argument("--output",     dest="output",  action="store_true",    help="enable output.")
    parser.add_argument("--no-output",  dest="output",  action="store_false",   help="disable output.")
    parser.add_argument("--Asimov",     dest="asimov",  action="store_true",    help="enable asimov dataset.")
    parser.add_argument("--no-Asimov",  dest="asimov",  action="store_false",   help="disable asimov dataset.")
    parser.add_argument("--scanning",   dest="scanning", action="store_true",   help="enable scanning.")
    parser.add_argument("--no-scanning",dest="scanning", action="store_false",  help="disable scanning.")
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
    output      = args.output
    fileNo      = args.fileNo
    fitDim      = args.fitDim
    nuMass      = args.nuMass
    asimov      = args.asimov
    scanning    = args.scanning
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

    scale = target * 10.0 * 10.0 / (dist * dist)  ## due to the PDFs are generated based on 10 kpc CCSNe and 20kton LS...
    bkgscale = target

    ##### Information Output for Check #####
    logging.debug("\n========== Configuration of Fitter ==========\n")
    logging.debug(f"=== This is {fitDim}-dimensional fitting for absolute neutrino mass. ===\n")
    logging.debug(f"=== Supernova model: {model} {modelNo} \n")
    logging.debug(f"=== Supernova distance: {dist} kpc\n")
    logging.debug(f"=== Mass Ordering of dataset: {MO} \n")
    if asimov:
        logging.debug(f"=== This is a Asimov dataset fit \n")
    logging.debug(f"=== Channel: eES={eES}, IBD={IBD}, pES={pES}\n")
    logging.debug(f"=== Detector target mass: {target * 20} kton\n")
    logging.debug("\n=============================================\n")
    ########################################

    print(f"***** Current signal statistical scale factor is {scale}.\n")
    print(f"***** Current background statistical scale factor is {bkgscale}.\n")

    glow, ghig = scanner.load_C14()
    c14rate = 0
    if C14:
        if C14level == "low":
            c14rate = glow.Eval(Ethr)
        else:
            c14rate = ghig.Eval(Ethr)
        print(f"***** Carbon-14 rate with detection threshold {Ethr:.2f} MeV = {c14rate} /s.\n") 

    if Ethr > 0.2:
        Eibd = Eth
    else:
        Eibd = 0.2

    channels = {}
    if pES:
        channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fitEmax=5, fileNo=fileNo, dist=dist, exp=exp)
        channels["pES"].setC14rate(c14rate)
    if IBD:
        channels["IBD"] = channel("IBD", MO, model, modelNo, Eibd, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
        channels["IBD"].setC14rate(0)
    if eES:
        channels["eES"] = channel("eES", MO, model, modelNo, Eibd, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
        channels["eES"].setC14rate(0)

    for cha in channels.values():

        cha.setNevtPerFile(1e5)
        cha.setStartEvtId(startevt)
        cha.setEndEvtId(endevt)

        cha.setScale(scale)
        cha.setBkgScale(bkgscale)

        # 1D PDF setting:

        pdffile1dpath = os.getenv("PDF1DFILEPATHNEW")
        #cha.setDataPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_{cha.MH}_10kpc_{cha.name}_nuMass0.0eV_TEobs2dPDFintegral_JUNO.root")
        cha.setDataPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_PDF_{cha.name}_{cha.MH}_10kpc_Ethr{cha.Ethr:.2f}MeV_nuMass0.0eV_tmin-0.030tmax0.170s_integral_newXS_JUNO.root")
        cha._load_datapdf1D()
        #cha.setNOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDFintegral_JUNO.root")
        cha.setNOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_PDF_{cha.name}_NO_10kpc_Ethr{cha.Ethr:.2f}MeV_nuMass{nuMass:.1f}eV_tmin-0.030tmax0.170s_integral_newXS_JUNO.root")
        #cha.setIOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDFintegral_JUNO.root")
        cha.setIOPdfFile1D(f"{pdffile1dpath}{model}{modelNo}_PDF_{cha.name}_IO_10kpc_Ethr{cha.Ethr:.2f}MeV_nuMass{nuMass:.1f}eV_tmin-0.030tmax0.170s_integral_newXS_JUNO.root")
        cha._load_pdf1D()

        #N_1DNO = stat_tool.get_statROI_from1Dpdf(cha.pdfNOfile1D)
        #N_1DIO = stat_tool.get_statROI_from1Dpdf(cha.pdfIOfile1D)
        #print(f"{cha.name} channel: 1D NO event number in ROI = {N_1DNO}, IO event number in ROI = {N_1DIO}.")

        if fitDim == 2:
            pdffile2dpath = os.getenv("PDF2DFILEPATHNEW")
            cha.setDataPdfFile2D(f"{pdffile2dpath}{model}{modelNo}_nuePDF_{cha.MH}_10kpc_{cha.name}_nuMass0.0eV_TEobs2dPDF_JUNO_newXS.root")
            cha._load_datapdf2D()
            cha.setNOPdfFile2D(f"{pdffile2dpath}{model}{modelNo}_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_JUNO_newXS.root")
            cha.setIOPdfFile2D(f"{pdffile2dpath}{model}{modelNo}_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_JUNO_newXS.root")
            cha._load_pdf2D()

        # Set Data Files
        if not asimov:
            if fitDim == 1:
                cha.setDataFile1D(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data1d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_T1D.root")
                cha._load_data1D()
            if fitDim == 2:
                datapath = os.getenv("DATA2DFILEPATH")
                cha.setDataFile2D(f"{datapath}Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_JUNO.root")

                #if cha.name == "pES":
                #    if not C14:
                #        cha.setDataFile2D(f"{datapath}Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_JUNO_scale25.0.root")
                #        #cha.setDataFile2D(f"{datapath}Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_noC14.root")
                #    else:
                #        if C14level == "low":
                #            cha.setDataFile2D(f"{datapath}Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_C14low.root")
                #        elif C14level == "high":
                #            cha.setDataFile2D(f"{datapath}Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_C14high.root")
                #else:
                #    cha.setDataFile2D(f"{datapath}Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_JUNO.root")

                cha._load_data2D()    # could get 1D or 2D dataset from the fitting requirement.

        FITTING_EVENT_NUM =  cha.getNevtPerFile() # the sample number to run...
    

    if scanning and not asimov:

        dt_scan = np.arange(-0.01, 0.04, 0.001)
        #dt_scan = np.array([0.0])
        chi2_scan = np.zeros(len(dt_scan))

        for i, dt in enumerate(dt_scan):
            scan = 0
            for cha in channels.values():
                if MO == "NO":
                    scan += cha.calc_Asimov_NLL_NO2D(dt)
                elif MO == "IO":
                    scan += cha.calc_Asimov_NLL_IO2D(dt)
            chi2_scan[i] = 2 * scan

        arr = np.vstack((dt_scan, chi2_scan))

        target = 20 * target
        filename = f"/junofs/users/miaoyu/supernova/simulation/toyMC/scanRes/{model}{modelNo}_{dist}kpc_{target}kton_{MO}_pES{pES}eES{eES}IBD{IBD}_nuMass{nuMass:.1f}eV_{Ethr:.2f}MeV_fitTmin{fitTmin:.3f}sfitTmax{fitTmax:.3f}s_asimov2Dscan_{exp}_absMass.csv"
        np.savetxt(filename, arr.T)
        sys.exit(-1)

    if asimov:

        # asimov dataset test
        logging.debug("\n ========= FITTING ASIMOV DATASET ========= \n")
    

        if doFit:

            dt_arr, nll_arr, TbestFit, locMinFit, _, _, _ = scanner.scanning_asimov_chain_absoluteMass(channels.values(), MO, fitDim)
            if output:
                X = np.vstack((nuMass*np.ones_like(dt_arr), dt_arr, nll_arr))
                outfilename = f"{model}{modelNo}_{dist}kpc_fit{fitDim}D_midNLL_{MO}_fitMass{nuMass:.1f}eV_IBD{IBD}eES{eES}pES{pES}_NonlRes.csv"
                print(f"***** Output filename = {outfilename}.")
                with open(outfilename, "ab") as f:
                    np.savetxt(f, X.T)
        else:
            for nuMass in np.arange(0.0, 2.1, 0.1):
                for cha in channels.values():
                    cha.setNOPdfFile2D(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                    cha.setIOPdfFile2D(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                    cha._load_pdf2D()

                tmpnll = scanner.scanning_asimov2D_allchannels_MP([0], channels.values(), MO)[0]
                print(nuMass, tmpnll)

    
    elif (not asimov) and (not scanning):
        ## fit toyMC fitting
        if doFit:
            ## Do fitting
            if endevt == 0:
                SUB_EVENT_NUM = int(FITTING_EVENT_NUM)
            else:
                SUB_EVENT_NUM = int(endevt - startevt)
            Tbest, locMin = np.zeros(SUB_EVENT_NUM), np.zeros(SUB_EVENT_NUM)

            if endevt == 0:
                endevt = FITTING_EVENT_NUM


            for ievt in tqdm(range(startevt, endevt, 1)):

                dt_arr, nll_arr, TbestFit, locMinFit, _, _, _ = scanner.scanning_toyMC_chain_absoluteMass_2Dnew(channels.values(), MO,  ievt)
                Tbest[ievt-startevt] = TbestFit
                locMin[ievt-startevt] = locMinFit

                
            if output:
                df = pd.DataFrame({
                    "Tbest" : Tbest,
                    "locMin" : locMin,
                })

                target = target * 20
                if not C14:
                    C14label = "noC14"
                else:
                    if C14level == "low":
                        C14label = "lowC14"
                    elif C14level == "high":
                        C14label = "highC14"
                outfilename = f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}{modelNo}_{dist}kpc_{target}kton_{MO}_pES{pES}eES{eES}IBD{IBD}_nuMass{nuMass:.1f}eV_{Ethr:.2f}MeV_fitTmin{fitTmin:.3f}sfitTmax{fitTmax:.3f}s_{C14label}_start{startevt}end{endevt}_PoisToyDataTobs{fitDim:d}D_{exp}_absMass_limitdTtoAvoidZeroProb.csv"
                df.to_csv(outfilename)
                print(f"\n *********** Data written in {outfilename}")


if __name__ == "__main__":
    fit_absolute_mass()


