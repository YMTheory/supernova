import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import argparse
from tqdm import tqdm
from channel_analyser import channel

import ROOT

def scanning(dT_arr, channels, ievt, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            data = cha.get_one_event(ievt)
            if MO == 'NO':
                val += cha.calc_binnedNLL_NO(data, dT)
                #val += cha.calc_NLL_NO(data, dT)
            else:
                val += cha.calc_binnedNLL_IO(data, dT)
                #val += cha.calc_NLL_IO(data, dT)
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





if __name__ == "__main__":

    MO      = "NO"
    model   = "Garching82503"
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
    exp     = "JUNO"
    startevt = 0
    endevt  = 0
    asimov  = False
    C14bkg  = False
    C14level = "high"
    doFit = True
    
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
    dist    = args.dist
    exp     = args.exp
    startevt = args.start
    endevt  = args.end
    asimov  = args.asimov
    C14bkg  = args.C14bkg
    C14level = args.C14level
    doFit   = args.fit
    

    if C14bkg:
        glow, ghig = load_C14()

    channels = {}
    if pES:
        channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    if IBD:
        channels["IBD"] = channel("IBD", MO, model, modelNo, 0.20, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    if eES:
        channels["eES"] = channel("eES", MO, model, modelNo, 0.20, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    

    # Initialization, data loading...
    for cha in channels.values():
        cha.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/scale10/{model}{modelNo}_{cha.name}_data_{MO}_10kpc_thr{Ethr}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_merger.root")
        if not asimov:
            cha.setNevtPerFile(100000)
            cha.setStartEvtId(startevt)
            cha.setEndEvtId(endevt)
            cha._load_data()
        cha._load_pdf()

        FITTING_EVENT_NUM =  cha.getNevtPerFile() # the sample number to run...
    
    dT_arr  = np.arange(-10, 11, 1) # coarse scanning, the step size is the bin width of PDFs (1ms/step)
    
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
            else:
                bestNLL = locMinFitIO

            fig, ax = plt.subplots(figsize=(6, 4))
            x1 = np.arange(TbestFitNO-3, TbestFitNO+3, 0.01)
            ax.plot(x1, 2*(aNO*x1**2+bNO*x1+cNO) - 2*bestNLL, lw=2, color="blue", label="NO")
            x2 = np.arange(TbestFitIO-3, TbestFitIO+3, 0.01)
            ax.plot(x2, 2*(aIO*x2**2+bIO*x2+cIO) - 2*bestNLL, lw=2, color="red",  label="IO")

            if MO == "NO":
                ax.hlines(locMinFitIO*2-bestNLL*2, x2[0]+1, x2[-1]-1, linestyle="--", lw=2, color="black")
                ax.text(TbestFitIO, locMinFitIO*2-bestNLL*2-2, r"$\Delta \chi^2 =$"+f"{locMinFitIO*2-bestNLL*2:.1f}", fontsize=16)
            if MO == "IO":
                ax.hlines(locMinFitNO*2-bestNLL*2, x1[0]+1, x1[-1]-1, linestyle="--", lw=2, color="black")
                ax.text(TbestFitNO+0.5, locMinFitNO*2-bestNLL*2-2, r"$\Delta \chi^2 =$"+f"{locMinFitNO*2-bestNLL*2:.1f}", fontsize=16)

            ax.set_xlabel(r"$\Delta t$ [ms]", fontsize=19)
            ax.set_ylabel(r"$\chi^2$", fontsize=19)
            ax.set_ylim(0, 15)
            ax.legend(prop={"size":16}, frameon=True)
            ax.tick_params(axis="both", labelsize=18)
            plt.tight_layout()
            plt.savefig(f"./plots/{model}{modelNo}_data{MO}_{dist}kpc_AsimovFit.pdf")
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
                    nllNO_oneEvt = scanning(dT_arr, channels.values(), ievt, "NO")          # coarse scanning
                Tbest, locMin, _ = find_locMin(dT_arr, nllNO_oneEvt)
                dT_arr_fine = generate_fine_dTarr(Tbest)
                if C14bkg:
                    nllNO_oneEvt = scanning_withbkg(dT_arr_fine, channels.values(), ievt, "NO", C14level, Tbkg)
                else:
                    nllNO_oneEvt = scanning(dT_arr_fine, channels.values(), ievt, "NO")     # fine scanning
                TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllNO_oneEvt)

                TbestNO[ievt-startevt] = TbestFit
                locMinNO[ievt-startevt] = locMinFit
                #print("NO", TbestFit, locMinFit)
        
                
                ## Fitting w/ IO PDF
                if C14bkg:
                    nllIO_oneEvt = scanning_withbkg(dT_arr, channels.values(), ievt, "IO", C14level, Tbkg)
                else:
                    nllIO_oneEvt = scanning(dT_arr, channels.values(), ievt, "IO")          # coarse scanning
                Tbest, locMin, _ = find_locMin(dT_arr, nllIO_oneEvt)
                dT_arr_fine = generate_fine_dTarr(Tbest)
                if C14bkg:
                    nllIO_oneEvt = scanning_withbkg(dT_arr_fine, channels.values(), ievt, "IO", C14level, Tbkg)     # fine scanning
                else:
                    nllIO_oneEvt = scanning(dT_arr_fine, channels.values(), ievt, "IO")     # fine scanning
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
            df.to_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}{modelNo}_{dist}kpc_{MO}_{cha}_{Ethr:.2f}MeV_fitTmax{fitTmax:d}ms_fileNo{fileNo:d}_new_{exp}_start{startevt}end{endevt}_noFit_binning_scale10.csv")
        
