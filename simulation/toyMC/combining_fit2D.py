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



if __name__ == "__main__" :
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
    """

    logging.basicConfig(filename='example.log',level=logging.DEBUG)

    print("==================================\n")
    print(f"Total CPU number is {cpu_count()}")
    print("==================================\n")

    MO      = "NO"
    model   = "Garching"
    modelNo = 82703
    Ethr    = 0.10
    fitTmin = -0.02
    fitTmax = 0.02
    nuMass  = 0.0
    fileNo  = 0
    dist    = 10.
    target  = 1.
    exp     = "JUNO"
    startevt= 0
    endevt  = 100
    doFit   = True
    fitDim  = 2
    asimov  = False
    test    = True
    eES     = True
    IBD     = True
    pES     = True
    useMass = False
    C14     = False
    C14level="low"

    parser = argparse.ArgumentParser(description='Arguments of SNNu analyser.')
    parser.add_argument('--model',      type=str,       default="Garching",     help="Model name (option: Garching, Burrows).")
    parser.add_argument('--modelNo',    type=int,       default=82703,          help="Model number.")
    parser.add_argument('--MO',         type=str,       default="NO",           help="Mass ordering for the dataset.")
    parser.add_argument("--start",      type=int,       default=0,              help="Start event for fit.")
    parser.add_argument("--end",        type=int,       default=0,              help="End event for fit.")
    parser.add_argument("--fileNo",     type=int,       default=0,              help="File No.")
    parser.add_argument("--doFit",      dest="fit",     action="store_true",    help="enable fitting.")
    parser.add_argument("--no-doFit",   dest="fit",     action="store_false",   help="disable fitting.")
    parser.add_argument("--useMass",    dest="useMass", action="store_true",    help="enable mass effect.")
    parser.add_argument("--no-useMass", dest="useMass", action="store_false",   help="disable mass effect.")
    parser.add_argument("--Asimov",     dest="asimov",  action="store_true",    help="enable asimov dataset.")
    parser.add_argument("--no-Asimov",  dest="asimov",  action="store_false",   help="disable asimov dataset.")
    parser.add_argument("--test",       dest="test",    action="store_true",    help="enable asimov test.")
    parser.add_argument("--no-test",    dest="test",    action="store_false",   help="disable asimov test.")
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
    test        = args.test
    eES         = args.eES
    IBD         = args.IBD
    pES         = args.pES
    Ethr        = args.Ethr
    useMass     = args.useMass
    fitTmin     = args.fitTmin
    fitTmax     = args.fitTmax
    C14         = args.C14
    C14level    = args.C14level
    dist        = args.dist
    target      = args.target

    scale = target * 10.0 * 10.0 / (dist * dist)  ## due to the PDFs are generated based on 10 kpc CCSNe and 20kton LS...
    bkgscale = target

    if not useMass :
        nuMass = 0.0

    if not useMass:
        fitMode = "mass ordering"
    else:
        fitMode = "absolute mass"

    ##### Information Output for Check #####
    logging.debug("\n========== Configuration of Fitter ==========\n")
    logging.debug(f"=== This is {fitDim}-dimensional fitting for {fitMode} ===\n")
    logging.debug(f"=== Supernova model: {model} {modelNo} \n")
    logging.debug(f"=== Supernova distance: {dist} kpc\n")
    logging.debug(f"=== Mass Ordering of dataset: {MO} \n")
    if useMass:
        logging.debug(f"=== Fitting neutrino mass: {nuMass:.1f} eV \n")
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
        Eibd = Ethr
    else:
        Eibd = 0.2

    channels = {}
    if pES:
        channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fitEmax=5, fileNo=fileNo, dist=dist, exp=exp)
        channels["pES"].setC14rate(c14rate)
        if C14:
            channels["pES"].setNOPdfwithBkgFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_pES_nuMass0.0eV_TEobs2dPDF_rebin_c14{C14level}_new.root")
            channels["pES"].setIOPdfwithBkgFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_pES_nuMass0.0eV_TEobs2dPDF_rebin_c14{C14level}_new.root")
            channels["pES"]._load_pdf2DwithBkg()
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

        if useMass: # absolute mass fitting
            cha.setNOPdfFile0Path(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{cha.name}_NO_10kpc_Ethr{cha.Ethr:.2f}MeV_nuMass0.0_Tobs1D.root")  # -> only eES and IBD avaiable now
            cha.setIOPdfFile0Path(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{cha.name}_IO_10kpc_Ethr{cha.Ethr:.2f}MeV_nuMass0.0_Tobs1D.root")
            cha._load_pdf0()
            cha.setNOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{cha.name}_NO_10kpc_Ethr{cha.Ethr:.2f}MeV_nuMass{nuMass:.1f}_Tobs1D.root")
            cha.setIOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{cha.name}_IO_10kpc_Ethr{cha.Ethr:.2f}MeV_nuMass{nuMass:.1f}_Tobs1D.root")
            cha._load_pdf()
            if fitDim == 2:
                cha.setNOPdf2DFile0Path(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass0.0eV_TEobs2dPDF_v2.root")
                cha.setIOPdf2DFile0Path(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass0.0eV_TEobs2dPDF_v2.root")
                cha._load_pdf2D0()
                cha.setNOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                cha.setIOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                cha._load_data2D()
            
        else:   # MO fitting
            #cha.setNOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_NO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
            #cha.setIOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_IO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
            #cha.setNOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_NO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_THEIA100.root")
            #cha.setIOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_IO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_THEIA100.root")
            #cha.setNOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_NO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_HyperK.root")
            #cha.setIOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_IO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_HyperK.root")
            cha.setNOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_NO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
            cha.setIOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_IO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
            cha._load_pdf()
            if fitDim == 2:
                #cha.setNOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                #cha.setIOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                if cha.name == "pES":
                    if not C14:
                        cha.setNOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_JUNO_rebin.root")
                        cha.setIOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_JUNO_rebin.root")
                    else:
                        if C14level == "low":
                            cha.setNOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_rebin_c14low_new.root")
                            cha.setIOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_rebin_c14low_new.root")
                        elif C14level == "high":
                            cha.setNOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_rebin_c14high_new.root")
                            cha.setIOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_rebin_c14high_new.root")
                    cha._load_pdf2D()

            # Set Data Files
            if not asimov:
                cha.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data1d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_T1D.root")
                cha._load_data_ak()
                if fitDim == 2 and cha == "pES":
                    if not C14:
                        cha.setData2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_noC14.root")
                    else:
                        if C14level == "low":
                            cha.setData2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_C14low.root")
                        elif C14level == "high":
                            cha.setData2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_unbinnedData_{cha.MH}_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_TEobs2D_rebin_C14high.root")

                    cha._load_data2D()    # could get 1D or 2D dataset from the fitting requirement.

            FITTING_EVENT_NUM =  cha.getNevtPerFile() # the sample number to run...
    
    if test:
        logging.debug("\n ===== Run Asimov Test ===== \n")
        dT_arr = [0.0019999999999999896]
        if fitDim == 2:
            nllNO_coarse = scanner.scanning_asimov2D_withBkg_MP(dT_arr, channels.values(), "NO", "MO")
            nllIO_coarse = scanner.scanning_asimov2D_withBkg_MP(dT_arr, channels.values(), "IO", "MO")
        elif fitDim == 1:
            nllNO_coarse = scanner.scanning_asimov1D(dT_arr, channels.values(), "NO", "MO", C14, "low")
            nllIO_coarse = scanner.scanning_asimov1D(dT_arr, channels.values(), "IO", "MO", C14, "low")

        print(nllNO_coarse)
        print(nllIO_coarse)

    
    elif (not test) and asimov:
        # asimov dataset test
        logging.debug("\n ========= FITTING ASIMOV DATASET ========= \n")
        if not useMass:
            dT_arr = np.arange(-0.01, 0.011, 0.001)  # corase scanning, 1 ms scanning step
            #dT_arr = np.arange(-0.001, 0.001, 0.001)  # corase scanning, 1 ms scanning step
            if fitDim == 2:
                if C14:
                    nllNO_coarse = scanner.scanning_asimov2D_withBkg_MP(dT_arr, channels.values(), "NO", "MO")
                    #nllNO_coarse = scanning_asimov2D_withBkg(dT_arr, channels.values(), "NO", "MO")
                else:
                    nllNO_coarse = scanner.scanning_asimov2D_MP(dT_arr, channels.values(), "NO", "MO")
                    #nllNO_coarse = scanning_asimov2D(dT_arr, channels.values(), "NO", "MO")
                print("nllNO_coarse", nllNO_coarse)

            elif fitDim == 1:
                nllNO_coarse = scanner.scanning_asimov1D(dT_arr, channels.values(), "NO", "MO", C14, "low")
                #print(nllNO_coarse)
            Tbest, locMin, _ = scanner.find_locMin(dT_arr, nllNO_coarse)
            print("Rough scanning Tbest ", Tbest, "locMin ", locMin)
            dT_arr_fine = scanner.generate_fine_dTarr(Tbest)
            print(dT_arr_fine)
            if fitDim == 2:
                if C14:
                    nllNO_fine = scanner.scanning_asimov2D_withBkg_MP(dT_arr_fine, channels.values(), "NO", "MO")
                    #nllNO_fine = scanning_asimov2D_withBkg(dT_arr_fine, channels.values(), "NO", "MO")
                else:
                    nllNO_fine = scanner.scanning_asimov2D_MP(dT_arr_fine, channels.values(), "NO", "MO")
                    #nllNO_fine = scanning_asimov2D(dT_arr_fine, channels.values(), "NO", "MO")
                print("nllNO_fine", nllNO_fine)
            elif fitDim == 1:
                nllNO_fine = scanner.scanning_asimov1D(dT_arr_fine, channels.values(), "NO", "MO", C14, "low")     # fine scanning
            TbestFitNO, locMinFitNO, aNO, bNO, cNO = scanner.parabola_fit(dT_arr_fine, nllNO_fine, param=True)
            print(f"NO pdf fit {MO} Asimov data -> {TbestFitNO*1000} ms, {locMinFitNO}")

            if fitDim == 2:
                if C14:
                    nllIO_coarse = scanner.scanning_asimov2D_withBkg_MP(dT_arr, channels.values(), "IO", "MO")
                    #nllIO_coarse = scanning_asimov2D_withBkg(dT_arr, channels.values(), "IO", "MO")
                else:
                    nllIO_coarse = scanner.scanning_asimov2D_MP(dT_arr, channels.values(), "IO", "MO")
                    #nllIO_coarse = scanning_asimov2D(dT_arr, channels.values(), "IO", "MO")
                print("nllIO_coarse", nllIO_coarse)
            elif fitDim == 1:
                nllIO_coarse = scanner.scanning_asimov1D(dT_arr, channels.values(), "IO", "MO", C14, "low")
            #print(nllIO_coarse)
            Tbest, locMin, _ = scanner.find_locMin(dT_arr, nllIO_coarse)
            print("Rough scanning Tbest ", Tbest, "locMin ", locMin)
            dT_arr_fine = scanner.generate_fine_dTarr(Tbest)
            print(dT_arr_fine)
            if fitDim == 2:
                if C14:
                    nllIO_fine = scanner.scanning_asimov2D_withBkg_MP(dT_arr_fine, channels.values(), "IO", "MO")
                    #nllIO_fine = scanning_asimov2D_withBkg(dT_arr_fine, channels.values(), "IO", "MO")
                else:
                    nllIO_fine = scanner.scanning_asimov2D_MP(dT_arr_fine, channels.values(), "IO", "MO")
                    #nllIO_fine = scanning_asimov2D(dT_arr_fine, channels.values(), "IO", "MO")
                print("nllIO_fine", nllIO_fine)
            elif fitDim == 1:
                nllIO_fine = scanner.scanning_asimov1D(dT_arr_fine, channels.values(), "IO", "MO", C14, "low")     # fine scanning
            TbestFitIO, locMinFitIO, aIO, bIO, cIO = scanner.parabola_fit(dT_arr_fine, nllIO_fine, param=True)
            print(f"IO pdf fit {MO} Asimov data -> {TbestFitIO*1000} ms, {locMinFitIO}")

            print(f"Output {Ethr} {TbestFitNO*1000} {locMinFitNO} {TbestFitIO*1000} {locMinFitIO}")

            if MO == "NO":
                bestNLL = locMinFitNO
                #print(f"Asimov dataset sensitivity of NO data = {2*(locMinFitIO - bestNLL)}")
            else:
                bestNLL = locMinFitIO
                #print(f"Asimov dataset sensitivity of IO data = {2*(locMinFitNO - bestNLL)}")

        else:
            print(f"Fit {MO} Asimov dataset with different neutrino mass pdfs inputs.")
            dT_arr = np.arange(-0.01, 0.011, 0.001)
            if fitDim == 2:
                nllNO_coarse = scanner.scanning_asimov2D(dT_arr, channels.values(), MO, "Mass")
            elif fitDim == 1:
                nllNO_coarse = scanner.scanning_asimov1D(dT_arr, channels.values(), MO, "Mass")
            print(nllNO_coarse)
            Tbest, locMin, _ = scanner.find_locMin(dT_arr, nllNO_coarse)
            dT_arr_fine = scanner.generate_fine_dTarr(Tbest)
            if fitDim == 2:
                nllNO_fine = scanner.scanning_asimov2D(dT_arr_fine, channels.values(), MO, "Mass")     # fine scanning
            elif fitDim == 1:
                nllNO_fine = scanner.scanning_asimov1D(dT_arr_fine, channels.values(), MO, "Mass")     # fine scanning
            TbestFitNO, locMinFitNO, aNO, bNO, cNO = scanner.parabola_fit(dT_arr_fine, nllNO_fine, param=True)
            print(f"Absolute mass fit {MO} Asimov data -> {TbestFitNO} s, {locMinFitNO}")
            outfn = f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/Garching82703_10kpc_{MO}_eESonly_asimovFit{fitDim}D_nuMass{nuMass:.1f}eV.txt"
            with open(outfn, "w") as fo:
                fo.write(f"{TbestFitNO} {locMinFitNO}")


    else:
        if doFit:
            ## Do fitting
            dT_arr = np.arange(-0.01, 0.011, 0.001)  # corase scanning, 1 ms scanning step
            if endevt == 0:
                SUB_EVENT_NUM = int(FITTING_EVENT_NUM)
            else:
                SUB_EVENT_NUM = int(endevt - startevt)
            TbestNO, TbestIO = np.zeros(SUB_EVENT_NUM), np.zeros(SUB_EVENT_NUM)
            locMinNO, locMinIO = np.zeros(SUB_EVENT_NUM), np.zeros(SUB_EVENT_NUM)
            NeES = np.zeros(SUB_EVENT_NUM)
            NpES = np.zeros(SUB_EVENT_NUM)
            NIBD = np.zeros(SUB_EVENT_NUM)

            if endevt == 0:
                endevt = FITTING_EVENT_NUM


            for ievt in tqdm(range(startevt, endevt, 1)):
                if eES:
                    NeES[ievt-startevt] =  channels["eES"].getNsigCurrentEvent(ievt)
                if pES:
                    NpES[ievt-startevt] = channels["pES"].getNsigCurrentEvent(ievt)
                if IBD:
                    NIBD[ievt-startevt] = channels["IBD"].getNsigCurrentEvent(ievt)

                if fitDim == 2:
                    nllNO_oneEvt = scanner.scanning2D(dT_arr, channels.values(), ievt, "NO")
                elif fitDim == 1:
                    nllNO_oneEvt = scanner.scanning1D(dT_arr, channels.values(), ievt, "NO")
                ## print("Rough nllNO_oneEvt ", nllNO_oneEvt)
                Tbest, locMin, _ = scanner.find_locMin(dT_arr, nllNO_oneEvt)
                ## print(Tbest, locMin)
                dT_arr_fine = scanner.generate_fine_dTarr(Tbest)
                if fitDim == 2:
                    nllNO_oneEvt = scanner.scanning2D(dT_arr_fine, channels.values(), ievt, "NO")
                elif fitDim == 1:
                    nllNO_oneEvt = scanner.scanning1D(dT_arr_fine, channels.values(), ievt, "NO")
                ## print("Fine  nllNO_oneEvt ", nllNO_oneEvt)
                TbestFit, locMinFit = scanner.parabola_fit(dT_arr_fine, nllNO_oneEvt)
                ## print(TbestFit, locMinFit)
            
                TbestNO[ievt-startevt] = TbestFit
                locMinNO[ievt-startevt] = locMinFit
            
            for ievt in tqdm(range(startevt, endevt, 1)):
                if fitDim == 2:
                    nllIO_oneEvt = scanner.scanning2D(dT_arr, channels.values(), ievt, "IO")
                elif fitDim == 1:
                    nllIO_oneEvt = scanner.scanning1D(dT_arr, channels.values(), ievt, "IO")
                Tbest, locMin, _ = scanner.find_locMin(dT_arr, nllIO_oneEvt)
                dT_arr_fine = scanner.generate_fine_dTarr(Tbest)
                if fitDim == 2:
                    nllIO_oneEvt = scanner.scanning2D(dT_arr_fine, channels.values(), ievt, "IO")
                elif fitDim == 1:
                    nllIO_oneEvt = scanner.scanning1D(dT_arr_fine, channels.values(), ievt, "IO")
                TbestFit, locMinFit = scanner.parabola_fit(dT_arr_fine, nllIO_oneEvt)
            
                TbestIO[ievt-startevt] = TbestFit
                locMinIO[ievt-startevt] = locMinFit

                
            if MO == "NO":
                sens = 2 * (locMinIO - locMinNO)
            else:
                sens = -2 * (locMinIO - locMinNO)

            df = pd.DataFrame({
                "locMinNO" : locMinNO,
                "locMinIO" : locMinIO,
                "sens" : sens,
                "TbestNO" : TbestNO,
                "TbestIO" : TbestIO,
                "NeES"    : NeES,
                "NpES"    : NpES,
                "NIBD"    : NIBD,
            })
            target = target * 20
            if not C14:
                C14label = "noC14"
            else:
                if C14level == "low":
                    C14label = "lowC14"
                elif C14level == "high":
                    C14label = "highC14"
            df.to_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}{modelNo}_{dist}kpc_{target}kton_{MO}_pESeESIBD_useMass{useMass}_nuMass{nuMass:.1f}eV_{Ethr:.2f}MeV_fitTmin{fitTmin:.3f}sfitTmax{fitTmax:.3f}s_{C14label}_start{startevt}end{endevt}_PoisToyDataTobs{fitDim:d}D_{exp}.csv")
