import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm
import ROOT
import warnings
warnings.filterwarnings("ignore")
import logging

from channel_analyser import channel


## Issues to note
## time units:
##### 2DPDF / 1DPDF_v2: s
##### 2D data: ms

def scanning1D(dT_arr, channels, ievt, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            #dataT = cha.get_one_event1D(ievt)
            dataT = cha.get_one_event(ievt)
            if MO == "NO":
                val += cha.calc_NLL_NO(dataT, dT)
            else:
                val += cha.calc_NLL_IO(dataT, dT)
        nll[idx] = val
    return nll


def scanning2D(dT_arr, channels, ievt, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            dataT, dataE = cha.get_one_event2D(ievt)
            if MO == "NO":
                val += cha.calc_NLL_NO2D(dataT, dataE, dT)
            else:
                val += cha.calc_NLL_IO2D(dataT, dataE, dT)
        nll[idx] = val
    return nll


def scanning_asimov1D(dT_arr, channels, MO, ty, useC14, level):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                val += cha.calc_Asimov_NLL_NO(dT, ty)
            else:
                val += cha.calc_Asimov_NLL_IO(dT, ty)
        nll[idx] = val

    return nll


def scanning_asimov2D_withBkg(dT_arr, channels, MO, ty):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                if cha.name == "pES":
                    val += cha.calc_Asimov_NLL_NO2D_withBkg(dT)
                else:
                    val += cha.calc_Asimov_NLL_NO(dT, ty)
            else:
                if cha.name == "pES":
                    val += cha.calc_Asimov_NLL_IO2D_withBkg(dT)
                else:
                    val += cha.calc_Asimov_NLL_IO(dT, ty)


        nll[idx] = val
    return nll


def scanning_asimov2D(dT_arr, channels, MO, ty):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                if cha.name == "pES":
                    val += cha.calc_Asimov_NLL_NO2D(dT, ty)
            else:
                val += cha.calc_Asimov_NLL_IO2D(dT, ty)

        nll[idx] = val

    return nll

    
def find_locMin(dT_arr, nll_arr):
    idx = np.argmin(nll_arr)
    Tbest = dT_arr[idx]
    locMin = nll_arr[idx]
    return Tbest, locMin, idx


def generate_fine_dTarr(Tbest):
    TMIN = -9
    if Tbest < TMIN:
        Tbest = TMIN
    HALF_FITTING_NUM = 20
    TSTEP = 0.0001 # s
    dT_arr = np.zeros(2 * HALF_FITTING_NUM + 1)
    for i in range(2*HALF_FITTING_NUM+1):
        dT_arr[i] = Tbest - HALF_FITTING_NUM * TSTEP + i * TSTEP
    return dT_arr


def parabola_fit(dT_arr, nll_arr, param=False):
    Tbest, locMin, idx = find_locMin(dT_arr, nll_arr)
    ## check if Tbest at the edge:
    N = len(dT_arr)
    if idx < 2 or idx > N-3: # not enough points to fit
        return Tbest, locMin
    else:
        a, b, c = np.polyfit(dT_arr[idx-2:idx+3], nll_arr[idx-2:idx+3], 2)
        Tbest = - b / 2 / a
        locMin = (4*a*c - b**2) / 4 / a
        if not param:
            return Tbest, locMin
        else:
            return Tbest, locMin, a, b, c


def load_C14():
    f = ROOT.TFile("/junofs/users/miaoyu/supernova/production/PDFs/backgrounds/C14/C14_rate.root", "read")
    glow = f.Get("c14_low")
    ghig = f.Get("c14_high")

    return glow, ghig



if __name__ == "__main__" :

    logging.basicConfig(filename='example.log',level=logging.DEBUG)

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

    glow, ghig = load_C14()
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
            pES._load_pdf2DwithBkg()
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
                cha.setNOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_NO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                cha.setIOPdf2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_nuePDF_IO_10kpc_{cha.name}_nuMass{nuMass:.1f}eV_TEobs2dPDF_v2.root")
                cha._load_pdf2D()

            # Set Data Files
            #cha.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data1d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_binning_newv2.root")
            #cha.setData2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_new2D.root")
            #cha.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data1d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_binning_newv2.root")
            #cha.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data1d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_T1D.root")
            cha.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data1d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_T1D_THEIA100.root")
            #cha._load_data_ak()
            cha.setData2DFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_{cha.name}_unbinneddata_{cha.MH}_10.0kpc_thr{cha.Ethr:.2f}MeV_Tmin-20msTmax20ms_Tobs2D.root")
            #cha._load_data2D()    # could get 1D or 2D dataset from the fitting requirement.
            #elif fitDim == 1:
            #    cha._load_data_ak()

            FITTING_EVENT_NUM =  cha.getNevtPerFile() # the sample number to run...
    
    if asimov:
        # asimov dataset test
        logging.debug("\n ========= FITTING ASIMOV DATASET ========= \n")
        if not useMass:
            dT_arr = np.arange(-0.01, 0.011, 0.001)  # corase scanning, 1 ms scanning step
            if fitDim == 2:
                if C14:
                    nllNO_coarse = scanning_asimov2D_withBkg(dT_arr, channels.values(), "NO", "MO")
                else:
                    nllNO_coarse = scanning_asimov2D(dT_arr, channels.values(), "NO", "MO")

            elif fitDim == 1:
                nllNO_coarse = scanning_asimov1D(dT_arr, channels.values(), "NO", "MO", C14, "low")
                #print(nllNO_coarse)
            Tbest, locMin, _ = find_locMin(dT_arr, nllNO_coarse)
            dT_arr_fine = generate_fine_dTarr(Tbest)
            if fitDim == 2:
                if C14:
                    nllNO_fine = scanning_asimov2D_withBkg(dT_arr, channels.values(), "NO", "MO")
                else:
                    nllNO_fine = scanning_asimov2D(dT_arr, channels.values(), "NO", "MO")
            elif fitDim == 1:
                nllNO_fine = scanning_asimov1D(dT_arr_fine, channels.values(), "NO", "MO", C14, "low")     # fine scanning
            TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dT_arr_fine, nllNO_fine, param=True)
            #print(f"NO pdf fit {MO} Asimov data -> {TbestFitNO*1000} ms, {locMinFitNO}")

            if fitDim == 2:
                if C14:
                    nllIO_coarse = scanning_asimov2D_withBkg(dT_arr, channels.values(), "IO", "MO")
                else:
                    nllIO_coarse = scanning_asimov2D(dT_arr, channels.values(), "IO", "MO")
            elif fitDim == 1:
                nllIO_coarse = scanning_asimov1D(dT_arr, channels.values(), "IO", "MO", C14, "low")
            #print(nllIO_coarse)
            Tbest, locMin, _ = find_locMin(dT_arr, nllIO_coarse)
            dT_arr_fine = generate_fine_dTarr(Tbest)
            if fitDim == 2:
                if C14:
                    nllIO_fine = scanning_asimov2D_withBkg(dT_arr, channels.values(), "IO", "MO")
                else:
                    nllIO_fine = scanning_asimov2D(dT_arr, channels.values(), "IO", "MO")
            elif fitDim == 1:
                nllIO_fine = scanning_asimov1D(dT_arr_fine, channels.values(), "IO", "MO", C14, "low")     # fine scanning
            TbestFitIO, locMinFitIO, aIO, bIO, cIO = parabola_fit(dT_arr_fine, nllIO_fine, param=True)
            #print(f"IO pdf fit {MO} Asimov data -> {TbestFitIO*1000} ms, {locMinFitIO}")

            print(f"Output {modelNo} {fitTmax} {TbestFitNO*1000} {locMinFitNO} {TbestFitIO} {locMinFitIO}")

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
                nllNO_coarse = scanning_asimov2D(dT_arr, channels.values(), MO, "Mass")
            elif fitDim == 1:
                nllNO_coarse = scanning_asimov1D(dT_arr, channels.values(), MO, "Mass")
            print(nllNO_coarse)
            Tbest, locMin, _ = find_locMin(dT_arr, nllNO_coarse)
            dT_arr_fine = generate_fine_dTarr(Tbest)
            if fitDim == 2:
                nllNO_fine = scanning_asimov2D(dT_arr_fine, channels.values(), MO, "Mass")     # fine scanning
            elif fitDim == 1:
                nllNO_fine = scanning_asimov1D(dT_arr_fine, channels.values(), MO, "Mass")     # fine scanning
            TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dT_arr_fine, nllNO_fine, param=True)
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
                    nllNO_oneEvt = scanning2D(dT_arr, channels.values(), ievt, "NO")
                elif fitDim == 1:
                    nllNO_oneEvt = scanning1D(dT_arr, channels.values(), ievt, "NO")
                ## print("Rough nllNO_oneEvt ", nllNO_oneEvt)
                Tbest, locMin, _ = find_locMin(dT_arr, nllNO_oneEvt)
                ## print(Tbest, locMin)
                dT_arr_fine = generate_fine_dTarr(Tbest)
                if fitDim == 2:
                    nllNO_oneEvt = scanning2D(dT_arr_fine, channels.values(), ievt, "NO")
                elif fitDim == 1:
                    nllNO_oneEvt = scanning1D(dT_arr_fine, channels.values(), ievt, "NO")
                ## print("Fine  nllNO_oneEvt ", nllNO_oneEvt)
                TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllNO_oneEvt)
                ## print(TbestFit, locMinFit)
            
                TbestNO[ievt-startevt] = TbestFit
                locMinNO[ievt-startevt] = locMinFit
            
            for ievt in tqdm(range(startevt, endevt, 1)):
                if fitDim == 2:
                    nllIO_oneEvt = scanning2D(dT_arr, channels.values(), ievt, "IO")
                elif fitDim == 1:
                    nllIO_oneEvt = scanning1D(dT_arr, channels.values(), ievt, "IO")
                Tbest, locMin, _ = find_locMin(dT_arr, nllIO_oneEvt)
                dT_arr_fine = generate_fine_dTarr(Tbest)
                if fitDim == 2:
                    nllIO_oneEvt = scanning2D(dT_arr_fine, channels.values(), ievt, "IO")
                elif fitDim == 1:
                    nllIO_oneEvt = scanning1D(dT_arr_fine, channels.values(), ievt, "IO")
                TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllIO_oneEvt)
            
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

            df.to_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}{modelNo}_{dist}kpc_{target}kton_{MO}_pESeESIBD_useMass{useMass}_nuMass{nuMass:.1f}eV_{Ethr:.2f}MeV_fitTmin{fitTmin:.3f}sfitTmax{fitTmax:.3f}s_data1D_start{startevt}end{endevt}_PoisToyDataTobs{fitDim:d}D_{exp}.csv")
