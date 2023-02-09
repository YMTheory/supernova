import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt

import argparse
from tqdm import tqdm
from channel_analyser import channel

import ROOT

def scanning_binned(dT_arr, channels, ievt, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            data = cha.get_one_binned_event(ievt)
            if MO == "NO":
                val += cha.calc_binnedNLL_NO(data, dT)
            else:
                val += cha.calc_binnedNLL_IO(data, dT)
        nll[idx] = val
    return nll


def scanning(dT_arr, channels, ievt, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            data = cha.get_one_event(ievt)
            if MO == 'NO':
                val += cha.calc_NLL_NO(data, dT)
            else:
                val += cha.calc_NLL_IO(data, dT)
        nll[idx] = val 

    return nll


def scanning_withbkg(dT_arr, channels, ievt, MO, level, Tbkg):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            data = cha.get_one_event(ievt)
            if MO == 'NO':
                if cha.name == "pES":
                    val += cha.calc_NLL_NO_withbkg(data, dT, Tbkg)
                else:
                    val += cha.calc_NLL_NO(data, dT)
            else:
                if cha.name == "pES":
                    val += cha.calc_NLL_NO_withbkg(data, dT, Tbkg)
                else:
                    val += cha.calc_NLL_IO(data, dT)
        nll[idx] = val 

    return nll


def scanning_asimov(dT_arr, channels, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                val += cha.calc_Asimov_NLL_NO(dT)
            else:
                val += cha.calc_Asimov_NLL_IO(dT)

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
    TSTEP = 0.1 # ms
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


def generate_C14(level, glow, ghig):
    if level == "low":
        c14rate = glow.Eval(Ethr)
    elif level == "high":
        c14rate = ghig.Eval(Ethr)
    else:
        print("Unknown level description!")

    nc14 = c14rate * (fitTmax - fitTmin)* 1e-3
    nc14 = np.random.poisson(nc14)
    if nc14 <= 0:
        return []
    Tc14 = np.random.uniform(fitTmin, fitTmax, size=nc14)
    print(f"C14-rate = {c14rate} Hz with event number in window = {nc14}.")
    return Tc14

def fitDrawer(arr, tmin, tmax, pdfx, pdfy, dataMO, fitMO, evtNo, dT, C14level):
    fig, ax = plt.subplots(figsize=(6, 4))
    #ax.hist(arr, bins=tmax-tmin, range=(tmin, tmax), color='blue', label="Data")
    conts, _ = np.histogram(arr, bins=tmax-tmin, range=(tmin, tmax))
    cents = np.arange(tmin+0.5, tmax+0.5, 1)
    ax.errorbar(cents[conts>0], conts[conts>0], yerr=np.sqrt(conts[conts>0]), fmt="o", ms=6, lw=2, color="blue", label="Data")
    ax.plot(pdfx, pdfy, lw=2, color="black", label="Fit")
    ax.set_xlabel("t [ms]", fontsize=15)
    ax.set_ylabel("counts/ms", fontsize=15)
    ax.legend(prop={"size":13})
    ax.tick_params(axis="both", labelsize=14)
    ax.set_title(r"$\Delta t$" + f" = {dT:.2f} ms", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"./plots/fitExample_data{dataMO}_fit{fitMO}_evtNo{evtNo}_C14{C14level}.pdf")
    plt.show()

def fitDrawer_withbkg(arr, arrbkg, tmin, tmax, pdfx, pdfy, dataMO, fitMO, evtNo, dT, C14level, c14rate):
    fig, ax = plt.subplots(figsize=(6, 4))
    #ax.hist(arr, bins=tmax-tmin, range=(tmin, tmax), color='blue', label="Data")
    conts, _ = np.histogram(arr, bins=tmax-tmin, range=(tmin, tmax))
    cents = np.arange(tmin+0.5, tmax+0.5, 1)
    ax.errorbar(cents[conts>0], conts[conts>0], yerr=np.sqrt(conts[conts>0]), fmt="^", ms=6, lw=2, color="blue", label="Data (signals)")
    conts, _ = np.histogram(arrbkg, bins=tmax-tmin, range=(tmin, tmax))
    cents = np.arange(tmin+0.5, tmax+0.5, 1)
    ax.errorbar(cents[conts>0], conts[conts>0], yerr=np.sqrt(conts[conts>0]), fmt="v", ms=6, lw=2, color="orange", label="Data (backgrounds)")
    allarr = []
    for i in arr:
        allarr.append(i)
    for i in arrbkg:
        allarr.append(i)
    conts, _ = np.histogram(allarr, bins=tmax-tmin, range=(tmin, tmax))
    cents = np.arange(tmin+0.5, tmax+0.5, 1)
    ax.errorbar(cents[conts>0], conts[conts>0], yerr=np.sqrt(conts[conts>0]), fmt="o", ms=6, lw=2, color="green", fillstyle="none", label="Data (all)")
    ax.plot(pdfx, pdfy, lw=2, color="black", label="Fit (all)")
    ax.plot(pdfx, c14rate/1000.*np.ones(len(pdfx)), lw=2, linestyle="--", color="gray", label="Fit (bkg)")
    ax.set_xlabel("t [ms]", fontsize=15)
    ax.set_ylabel("counts/ms", fontsize=15)
    ax.legend(prop={"size":13}, ncol=2, frameon=True)
    ax.tick_params(axis="both", labelsize=14)
    ax.set_title(r"$\Delta t$" + f" = {dT:.2f} ms" + f" data{dataMO}+pdf{fitMO}+evtNO{evtNo}", fontsize=15)
    plt.tight_layout()
    plt.savefig(f"./plots/fitExample_data{dataMO}_fit{fitMO}_evtNo{evtNo}_C14{C14level}.pdf")
    #plt.show()







if __name__ == "__main__":

    MO      = "NO"
    model   = "Garching"
    modelNo = 82503
    Ethr    = 0.15
    output  = True
    fitTmin = 10
    fitTmax = 50
    fileNo  = 0
    dist    = 10
    eES     = True
    IBD     = True
    pES     = True
    CEvNS   = False
    exp     = "JUNO"
    startevt = 0
    endevt  = 0
    asimov  = False
    C14bkg  = False
    C14level = "high"
    doFit   = True
    binning = False
    binwidth = 0.01
    
    parser = argparse.ArgumentParser(description='Arguments of SNNu analyser.')
    parser.add_argument('--model', type=str, default='Garching', help='Model name of SNNu.')
    parser.add_argument('--modelNo', type=int, default=82503, help="modelNo")
    parser.add_argument('--MO', type=str, default="NO", help="Mass ordering for the dataset.")
    parser.add_argument('--Ethr' , type=float, default=0.20, help="Detection threshold for pES channel, unit MeV.")
    parser.add_argument('--output', dest='output', action="store_true", help="output csv file.")
    parser.add_argument('--no-output', dest='output',action="store_false", help="do not output csv file.")
    parser.add_argument('--fitTmin', type=int, default=10, help="Minimum fitting time.")
    parser.add_argument('--fitTmax', type=int, default=50, help="Maximum fitting time.")
    parser.add_argument('--fileNo', type=int, default=0, help="data file number.")
    parser.add_argument('--dist', type=int, default=10, help="CCSNe distance")
    parser.add_argument('--eES', dest='eES', action="store_true", help="enable eES")
    parser.add_argument('--no-eES', dest='eES', action="store_false", help="disable eES")
    parser.add_argument('--IBD', dest='IBD', action="store_true", help="enable IBD")
    parser.add_argument('--no-IBD', dest='IBD', action="store_false", help="disable IBD")
    parser.add_argument('--pES', dest='pES', action="store_true", help="enable pES")
    parser.add_argument('--no-pES', dest='pES', action="store_false", help="disable pES")
    parser.add_argument('--CEvNS', dest='CEvNS', action="store_true", help="enable CEvNS")
    parser.add_argument('--no-CEvNS', dest='CEvNS', action="store_false", help="disable CEvNS")
    parser.add_argument("--exp", type=str, default="", help="Experiment name")
    parser.add_argument("--start", type=int, default=0, help="Start event for fit.")
    parser.add_argument("--end", type=int, default=0, help="End event for fit.")
    parser.add_argument("--asimov", dest="asimov", action="store_true", help="enable Asimov dataset fit.")
    parser.add_argument("--no-asimov", dest="asimov", action="store_false", help="disable Asimov dataset fit.")
    parser.add_argument("--C14bkg", dest="C14bkg", action="store_true", help="enable C14 background.")
    parser.add_argument("--no-C14bkg", dest="C14bkg", action="store_false", help="disable C14 background.")
    parser.add_argument("--C14level", type=str, default="high", help="C14 level.")
    parser.add_argument("--doFit", dest="fit", action="store_true", help="enable time fitting.")
    parser.add_argument("--no-doFit", dest="fit", action="store_false", help="disable time fitting.")
    parser.add_argument("--binning", dest="binning", action="store_true", help="enable binning fitting.")
    parser.add_argument("--no-binning", dest="binning", action="store_false", help="disable binning fitting.")
    parser.add_argument('--binwidth' , type=float, default=0.01, help="Bin width for binning fit.")
    parser.add_argument("--drawFit", dest="drawFit", action="store_true", help="draw fit example.")
    parser.add_argument("--no-drawFit", dest="drawFit", action="store_false", help="not draw fit example.")
    parser.add_argument("--drawNo", type=int, default=0, help="event number for drawing fit example.")
    args = parser.parse_args()
    
    model   = args.model
    modelNo = args.modelNo
    Ethr    = args.Ethr
    MO      = args.MO
    output  = args.output
    fitTmin = args.fitTmin
    fitTmax = args.fitTmax
    fileNo  = args.fileNo
    eES     = args.eES
    IBD     = args.IBD
    pES     = args.pES
    CEvNS   = args.CEvNS
    dist    = args.dist
    exp     = args.exp
    startevt = args.start
    endevt  = args.end
    asimov  = args.asimov
    C14bkg  = args.C14bkg
    C14level = args.C14level
    doFit   = args.fit
    binning = args.binning
    binwidth = args.binwidth
    glow, ghig = None, None
    drawFit = args.drawFit
    drawNo  = args.drawNo
    

    if C14bkg:
        glow, ghig = load_C14()
        if C14level == "low":
            c14rate = glow.Eval(Ethr)
        elif C14level == "high":
            c14rate = ghig.Eval(Ethr)

    channels = {}
    if pES:
        channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
        if C14bkg:
            channels["pES"].c14rate = c14rate
    if IBD:
        channels["IBD"] = channel("IBD", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    if eES:
        channels["eES"] = channel("eES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    if CEvNS:
        channels["CEvNS"] = channel("CEvNS", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    

    # Initialization, data loading...
    for cha in channels.values():
        if not asimov:
            if drawFit:
                cha.setNevtPerFile(1e5)
                cha.setStartEvtId(drawNo)
                cha.setEndEvtId(drawNo+1)
                if not binning:
                    cha.setDataFilePath(f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha.name}_binneddata_{MO}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning.root")
                    cha._load_data_ak()
                if binning:
                    cha.setDataFilePath(f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha.name}_binneddata_{MO}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning_new.root")
                    cha._load_data_binned()
            else:
                cha.setNevtPerFile(1e5)
                cha.setStartEvtId(startevt)
                cha.setEndEvtId(endevt)
                if not binning:
                    #cha.setDataFilePath(f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha.name}_binneddata_{MO}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning.root")
                    cha.setDataFilePath(f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha.name}_binneddata_{MO}_{dist:.1f}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning.root")
                    cha._load_data_ak()
                if binning:
                    cha.setDataFilePath(f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha.name}_binneddata_{MO}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning_new.root")
                    cha._load_data_binned()
        else:
            cha.binwidth = binwidth
        
        cha._load_pdf()

        FITTING_EVENT_NUM =  cha.getNevtPerFile() # the sample number to run...
    
    dT_arr  = np.arange(-10, 11, 1) # coarse scanning, the step size is the bin width of PDFs (1ms/step)
    
    if drawFit:
        if C14bkg:
            Tbkg = generate_C14(C14level, glow, ghig)
            nllNO_oneEvt = scanning_withbkg(dT_arr, channels.values(), drawNo, "NO", C14level, Tbkg)
        else:
            if not binning:
                nllNO_oneEvt = scanning(dT_arr, channels.values(), drawNo, "NO")
            else:
                nllNO_oneEvt = scanning_binned(dT_arr, channels.values(), drawNo, "NO")
        Tbest, locMin, _ = find_locMin(dT_arr, nllNO_oneEvt)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        if C14bkg:
            nllNO_oneEvt = scanning_withbkg(dT_arr_fine, channels.values(), drawNo, "NO", C14level, Tbkg)     # fine scanning
        else:
            if not binning:
                nllNO_oneEvt = scanning(dT_arr_fine, channels.values(), drawNo, "NO")     # fine scanning
            else:
                nllNO_oneEvt = scanning_binned(dT_arr_fine, channels.values(), drawNo, "NO")     # fine scanning
        TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllNO_oneEvt)
        
        if C14bkg:
            x = np.arange(fitTmin, fitTmax, 1)
            y = channels["pES"]._pdfNO_func(x+TbestFit)
            y += c14rate / 1000.
            fitevt = channels["pES"].get_one_event(drawNo)
            fitDrawer_withbkg(fitevt, Tbkg, fitTmin, fitTmax, x, y, MO, "NO", drawNo, TbestFit, C14level, c14rate)
                
        else:
            x = np.arange(fitTmin, fitTmax, 1)
            y = channels["pES"]._pdfNO_func(x+TbestFit)
            fitevt = channels["pES"].get_one_event(drawNo)
            fitDrawer(fitevt, fitTmin, fitTmax, x, y, MO, "NO", drawNo, TbestFit, "No")
        

        ## Fitting w/ IO PDF
        if C14bkg:
            nllIO_oneEvt = scanning_withbkg(dT_arr, channels.values(), drawNo, "IO", C14level, Tbkg)
        else:
            if not binning:
                nllIO_oneEvt = scanning(dT_arr, channels.values(), drawNo, "IO")          # coarse scanning
            else:
                nllIO_oneEvt = scanning_binned(dT_arr, channels.values(), drawNo, "IO")
        Tbest, locMin, _ = find_locMin(dT_arr, nllIO_oneEvt)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        if C14bkg:
            nllIO_oneEvt = scanning_withbkg(dT_arr_fine, channels.values(), drawNo, "IO", C14level, Tbkg)     # fine scanning
        else:
            if not binning:
                nllIO_oneEvt = scanning(dT_arr_fine, channels.values(), drawNo, "IO")     # fine scanning
            else:
                nllIO_oneEvt = scanning_binned(dT_arr_fine, channels.values(), drawNo, "IO")     # fine scanning
        TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllIO_oneEvt)
    
        if C14bkg:
            x = np.arange(fitTmin, fitTmax, 1)
            y = channels["pES"]._pdfIO_func(x+TbestFit)
            y += c14rate / 1000.
            fitevt = channels["pES"].get_one_event(drawNo)
            fitDrawer_withbkg(fitevt, Tbkg, fitTmin, fitTmax, x, y, MO, "IO", drawNo, TbestFit, C14level, c14rate)
                
        else:
            x = np.arange(fitTmin, fitTmax, 1)
            y = channels["pES"]._pdfIO_func(x+TbestFit)
            fitevt = channels["pES"].get_one_event(drawNo)
            fitDrawer(fitevt, fitTmin, fitTmax, x, y, MO, "IO", drawNo, TbestFit, "No")

        sys.exit(-1)
        
    if asimov == True:
        print("=====> Fitting Asimov dataset <=====")
        if doFit:
        
            nllNO_coarse = scanning_asimov(dT_arr, channels.values(), "NO")
            Tbest, locMin, _ = find_locMin(dT_arr, nllNO_coarse)
            dT_arr_fine = generate_fine_dTarr(Tbest)
            nllNO_fine = scanning_asimov(dT_arr_fine, channels.values(), "NO")     # fine scanning
            TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dT_arr_fine, nllNO_fine, param=True)
            print(f"NO pdf fit {MO} Asimov data -> {TbestFitNO}ms, {locMinFitNO}")

            nllIO_coarse = scanning_asimov(dT_arr, channels.values(), "IO")
            print(nllIO_coarse)
            Tbest, locMin, _ = find_locMin(dT_arr, nllIO_coarse)
            dT_arr_fine = generate_fine_dTarr(Tbest)
            nllIO_fine = scanning_asimov(dT_arr_fine, channels.values(), "IO")     # fine scanning
            TbestFitIO, locMinFitIO, aIO, bIO, cIO = parabola_fit(dT_arr_fine, nllIO_fine, param=True)
            print(f"IO pdf fit {MO} Asimov data -> {TbestFitIO}ms, {locMinFitIO}")

            if MO == "NO":
                bestNLL = locMinFitNO
                print(f"Asimov dataset sensitivity of NO data = {2*(locMinFitIO - bestNLL)}")
            else:
                bestNLL = locMinFitIO
                print(f"Asimov dataset sensitivity of IO data = {2*(locMinFitNO - bestNLL)}")


            if 1:
                print("Draw asimov dataset plots...")
                fig, ax = plt.subplots(figsize=(6, 4))
                x1 = np.arange(TbestFitNO-4, TbestFitNO+4, 0.01)
                ax.plot(x1, 2*(aNO*x1**2+bNO*x1+cNO) - 2*bestNLL, lw=2, color="blue", label="NO")
                #for i in x1:
                #    print("NO", i, 2*(aNO*i**2+bNO*i+cNO)-2*bestNLL)
                x2 = np.arange(TbestFitIO-3, TbestFitIO+3, 0.01)
                ax.plot(x2, 2*(aIO*x2**2+bIO*x2+cIO) - 2*bestNLL, lw=2, color="red",  label="IO")
                #for i in x2:
                #    print("IO", i, 2*(aIO*i**2+bIO*i+cIO)-2*bestNLL)

                if MO == "NO":
                    ax.hlines(locMinFitIO*2-bestNLL*2, x2[0]+1, x2[-1]-1, linestyle="--", lw=2, color="black")
                    ax.text(TbestFitIO, locMinFitIO*2-bestNLL*2-2, r"$\Delta \chi^2 =$"+f"{locMinFitIO*2-bestNLL*2:.1f}", fontsize=16)
                    ax.vlines(TbestFitIO, 0, locMinFitIO*2-2*bestNLL+2, linestyle=":", lw=2, color="gray")
                    ax.text(TbestFitIO+0.5, 2, f"{TbestFitIO:.2f} ms", fontsize=16)
                if MO == "IO":
                    ax.hlines(locMinFitNO*2-bestNLL*2, x1[0]+1, x1[-1]-1, linestyle="--", lw=2, color="black")
                    ax.text(TbestFitNO+0.5, locMinFitNO*2-bestNLL*2-2, r"$\Delta \chi^2 =$"+f"{locMinFitNO*2-bestNLL*2:.1f}", fontsize=16)
                    ax.vlines(TbestFitNO, 0, locMinFitNO*2-2*bestNLL+2, linestyle=":", lw=2, color="gray")
                    ax.text(TbestFitNO+0.5, 2, f"{TbestFitNO:.2f} ms", fontsize=16)

                ax.set_xlabel(r"$\Delta t$ [ms]", fontsize=19)
                ax.set_ylabel(r"$\chi^2$", fontsize=19)
                ax.set_ylim(0, 15)
                ax.legend(loc="upper center", prop={"size":16}, frameon=True, ncol=2)
                ax.tick_params(axis="both", labelsize=18)
                plt.tight_layout()
                plt.savefig(f"./plots/{model}{modelNo}_data{MO}_{dist}kpc_{Ethr:.2f}MeV_{binwidth:.2f}ms_AsimovFit.pdf")
                plt.show()

        if not doFit:
            nllMinNO_oneEvt = scanning_asimov([0], channels.values(), "NO")[0]
            nllMinIO_oneEvt = scanning_asimov([0], channels.values(), "IO")[0]
            if MO == "NO":
                sens = 2 * (nllMinIO_oneEvt - nllMinNO_oneEvt)
            else:
                sens = 2 * (nllMinNO_oneEvt - nllMinIO_oneEvt)
            print(f"{MO} Asimov data, sensitivity = {sens:.2f}")

    if asimov == False:

        if endevt == 0:
            SUB_EVENT_NUM = FITTING_EVENT_NUM
        else:
            SUB_EVENT_NUM = endevt - startevt
        TbestNO, TbestIO = np.zeros(SUB_EVENT_NUM), np.zeros(SUB_EVENT_NUM)
        locMinNO, locMinIO = np.zeros(SUB_EVENT_NUM), np.zeros(SUB_EVENT_NUM)

        
        ### LOOP in all FITTING DATA:
        if endevt == 0:
            endevt = FITTING_EVENT_NUM
        for ievt in tqdm(range(startevt, endevt, 1)):
        #for ievt in tqdm(range(FITTING_EVENT_NUM)):
            
            if doFit:
                ## Fitting w/ NO PDF
                if C14bkg:
                    Tbkg = generate_C14(C14level, glow, ghig)
                    nllNO_oneEvt = scanning_withbkg(dT_arr, channels.values(), ievt, "NO", C14level, Tbkg)
                else:
                    if not binning:
                        nllNO_oneEvt = scanning(dT_arr, channels.values(), ievt, "NO")          # coarse scanning
                    else:
                        nllNO_oneEvt = scanning_binned(dT_arr, channels.values(), ievt, "NO")
                Tbest, locMin, _ = find_locMin(dT_arr, nllNO_oneEvt)
                dT_arr_fine = generate_fine_dTarr(Tbest)
                if C14bkg:
                    nllNO_oneEvt = scanning_withbkg(dT_arr_fine, channels.values(), ievt, "NO", C14level, Tbkg)
                else:
                    if not binning:
                        nllNO_oneEvt = scanning(dT_arr_fine, channels.values(), ievt, "NO")     # fine scanning
                    else:
                        nllNO_oneEvt = scanning_binned(dT_arr_fine, channels.values(), ievt, "NO")
                TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllNO_oneEvt)

                TbestNO[ievt-startevt] = TbestFit
                locMinNO[ievt-startevt] = locMinFit
                #print("NO", TbestFit, locMinFit)
        
                
                ## Fitting w/ IO PDF
                if C14bkg:
                    nllIO_oneEvt = scanning_withbkg(dT_arr, channels.values(), ievt, "IO", C14level, Tbkg)
                else:
                    if not binning:
                        nllIO_oneEvt = scanning(dT_arr, channels.values(), ievt, "IO")          # coarse scanning
                    else:
                        nllIO_oneEvt = scanning_binned(dT_arr, channels.values(), ievt, "IO")
                Tbest, locMin, _ = find_locMin(dT_arr, nllIO_oneEvt)
                dT_arr_fine = generate_fine_dTarr(Tbest)
                if C14bkg:
                    nllIO_oneEvt = scanning_withbkg(dT_arr_fine, channels.values(), ievt, "IO", C14level, Tbkg)     # fine scanning
                else:
                    if not binning:
                        nllIO_oneEvt = scanning(dT_arr_fine, channels.values(), ievt, "IO")     # fine scanning
                    else:
                        nllIO_oneEvt = scanning_binned(dT_arr_fine, channels.values(), ievt, "IO")     # fine scanning
                TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllIO_oneEvt)

                TbestIO[ievt-startevt] = TbestFit
                locMinIO[ievt-startevt] = locMinFit
                #print("IO", TbestFit, locMinFit)

            if not doFit:
                # no time smearing
                nllMinNO_oneEvt = scanning([0], channels.values(), ievt, "NO")[0]
                nllMinIO_oneEvt = scanning([0], channels.values(), ievt, "IO")[0]
                locMinNO[ievt - startevt] = nllMinNO_oneEvt
                TbestNO[ievt - startevt] = 0
                locMinIO[ievt - startevt] = nllMinIO_oneEvt
                TbestIO[ievt - startevt] = 0
        


        if MO == "NO":
            sens = 2 * (locMinIO - locMinNO)
        else:
            sens = 2 * (locMinNO - locMinIO)

        if output:
            df = pd.DataFrame({ 
                                "locMinNO" : locMinNO,
                                "locMinIO" : locMinIO,
                                "sens" : sens,
                                "TbestNO" : TbestNO,
                                "TbestIO" : TbestIO
                             })

            cha = ""
            if pES:
                cha += "pES"
            if eES:
                cha += "eES"
            if IBD:
                cha += "IBD"
            if CEvNS:
                cha += "CEvNS"
            
            if binning:
                binway = "binned"
            else:
                binway = "unbinned"

            if not C14bkg:
                bkgflag = "noC14bkg"
            else:
                bkgflag = "C14bkg" + C14level


            df.to_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}{modelNo}_{dist}kpc_{MO}_{cha}_{Ethr:.2f}MeV_fitTmax{fitTmax:d}ms_start{startevt}end{endevt}_{doFit}_{binway}_{bkgflag}_PoisToyData.csv")
        
