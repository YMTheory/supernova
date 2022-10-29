import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import argparse
from tqdm import tqdm
from channel_analyser import channel


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


def find_locMin(dT_arr, nll_arr):
    idx = np.argmin(nll_arr)
    Tbest = dT_arr[idx]
    locMin = nll_arr[idx]
    return Tbest, locMin, idx


def generate_fine_dTarr(Tbest):
    TMIN = -9
    if Tbest < TMIN:
        Tbest = TMIN
    HALF_FITTING_NUM = 10
    TSTEP = 0.1 # ms
    dT_arr = np.zeros(2 * HALF_FITTING_NUM + 1)
    for i in range(2*HALF_FITTING_NUM+1):
        dT_arr[i] = Tbest - HALF_FITTING_NUM * TSTEP + i * TSTEP
    return dT_arr


def parabola_fit(dT_arr, nll_arr):
    Tbest, locMin, idx = find_locMin(dT_arr, nll_arr)
    ## check if Tbest at the edge:
    N = len(dT_arr)
    if idx < 2 or idx > N-3: # not enough points to fit
        return Tbest, locMin
    else:
        a, b, c = np.polyfit(dT_arr[idx-2:idx+3], nll_arr[idx-2:idx+3], 2)
        Tbest = - b / 2 / a
        locMin = (4*a*c - b**2) / 4 / a
        return Tbest, locMin
    

if __name__ == "__main__":

    MO      = "NO"
    model   = "Garching82503"
    modelNo = 82503
    Ethr    = 0.15
    output  = True
    fitTmin = 10
    fitTmax = 50
    fileNo  = 0
    eES     = True
    IBD     = True
    pES     = True
    
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
    parser.add_argument('--eES', dest='eES', action="store_true", help="enable eES")
    parser.add_argument('--no-eES', dest='eES', action="store_false", help="disable eES")
    parser.add_argument('--IBD', dest='IBD', action="store_true", help="enable IBD")
    parser.add_argument('--no-IBD', dest='IBD', action="store_false", help="disable IBD")
    parser.add_argument('--pES', dest='pES', action="store_true", help="enable pES")
    parser.add_argument('--no-pES', dest='pES', action="store_false", help="disable pES")
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
    
    channels = {}
    if pES:
        channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo)
    if IBD:
        channels["IBD"] = channel("IBD", MO, model, modelNo, 0.20, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo)
    if eES:
        channels["eES"] = channel("eES", MO, model, modelNo, 0.20, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo)
    

    # Initialization, data loading...
    for cha in channels.values():
        cha._load_data()
        cha._load_pdf()


    FITTING_EVENT_NUM =  1000 # the sample number to run...
    dT_arr  = np.arange(-10, 11, 1) # coarse scannig, the step size is the bin width of PDFs (1ms/step)
    
    TbestNO, TbestIO = np.zeros(FITTING_EVENT_NUM), np.zeros(FITTING_EVENT_NUM)
    locMinNO, locMinIO = np.zeros(FITTING_EVENT_NUM), np.zeros(FITTING_EVENT_NUM)
    
    ### LOOP in all FITTING DATA:
    for ievt in tqdm(range(FITTING_EVENT_NUM)):
        
        ## Fitting w/ NO PDF
        nllNO_oneEvt = scanning(dT_arr, channels.values(), ievt, "NO")          # coarse scanning
        Tbest, locMin, _ = find_locMin(dT_arr, nllNO_oneEvt)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        nllNO_oneEvt = scanning(dT_arr_fine, channels.values(), ievt, "NO")     # fine scanning
        TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllNO_oneEvt)

        TbestNO[ievt] = TbestFit
        locMinNO[ievt] = locMinFit
    
        
        ## Fitting w/ IO PDF
        nllIO_oneEvt = scanning(dT_arr, channels.values(), ievt, "IO")          # coarse scanning
        Tbest, locMin, _ = find_locMin(dT_arr, nllIO_oneEvt)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        nllIO_oneEvt = scanning(dT_arr_fine, channels.values(), ievt, "IO")     # fine scanning
        TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllIO_oneEvt)

        TbestIO[ievt] = TbestFit
        locMinIO[ievt] = locMinFit


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

        df.to_csv(f"/junofs/users/miaoyu/supernova/analysis/MH_new/results/{model}{modelNo}_10kpc_{MO}_pESeESIBD_{Ethr:.2f}MeV_fitTmax{fitTmax:d}ms_fileNo{fileNo:d}.csv")
        
