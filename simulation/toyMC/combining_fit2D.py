import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm
import ROOT

from channel_analyser import channel


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


if __name__ == "__main__" :

    MO = "NO"
    model = "Garching"
    modelNo = 82703
    Ethr = 0.10
    fitTmin = -20
    fitTmax = 30
    fileNo = 0
    dist = 10
    exp = "JUNO"
    startevt = 0
    endevt  = 100

    parser = argparse.ArgumentParser(description='Arguments of SNNu analyser.')
    parser.add_argument('--MO', type=str, default="NO", help="Mass ordering for the dataset.")
    parser.add_argument("--start", type=int, default=0, help="Start event for fit.")
    parser.add_argument("--end", type=int, default=0, help="End event for fit.")
    args = parser.parse_args()

    MO = args.MO
    startevt = args.start
    endevt = args.end
    
    channels = {}
    channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    channels["IBD"] = channel("IBD", MO, model, modelNo, 0.20, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    channels["eES"] = channel("eES", MO, model, modelNo, 0.20, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)

    for cha in channels.values():

        cha.setNevtPerFile(1e5)
        cha.setStartEvtId(startevt)
        cha.setEndEvtId(endevt)

        cha._load_data2D()
        cha._load_pdf()
        cha._load_pdf2D()
        
        FITTING_EVENT_NUM =  cha.getNevtPerFile() # the sample number to run...
    
    dT_arr  = np.arange(-10, 11, 1) # coarse scanning, the step size is the bin width of PDFs (1ms/step)
    if endevt == 0:
        SUB_EVENT_NUM = FITTING_EVENT_NUM
    else:
        SUB_EVENT_NUM = endevt - startevt
    TbestNO, TbestIO = np.zeros(SUB_EVENT_NUM), np.zeros(SUB_EVENT_NUM)
    locMinNO, locMinIO = np.zeros(SUB_EVENT_NUM), np.zeros(SUB_EVENT_NUM)
    if endevt == 0:
        endevt = FITTING_EVENT_NUM


    #print(channels["pES"].get_one_event2D(0))
    #print(channels["eES"].get_one_event2D(0))
    #print(channels["IBD"].get_one_event2D(0))

    for ievt in tqdm(range(startevt, endevt, 1)):
        nllNO_oneEvt = scanning2D(dT_arr, channels.values(), ievt, "NO")
        #print("Rough nllNO_oneEvt ", nllNO_oneEvt)
        Tbest, locMin, _ = find_locMin(dT_arr, nllNO_oneEvt)
        #print(Tbest, locMin)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        nllNO_oneEvt = scanning2D(dT_arr_fine, channels.values(), ievt, "NO")
        #print("Fine  nllNO_oneEvt ", nllNO_oneEvt)
        TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllNO_oneEvt)
        #print(TbestFit, locMinFit)
    
        TbestNO[ievt-startevt] = TbestFit
        locMinNO[ievt-startevt] = locMinFit
    
    for ievt in tqdm(range(startevt, endevt, 1)):
        nllIO_oneEvt = scanning2D(dT_arr, channels.values(), ievt, "IO")
        Tbest, locMin, _ = find_locMin(dT_arr, nllIO_oneEvt)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        nllIO_oneEvt = scanning2D(dT_arr_fine, channels.values(), ievt, "IO")
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
    })

    df.to_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}{modelNo}_{dist}kpc_{MO}_pESeESIBD_{Ethr:.2f}MeV_fitTmax{fitTmax:d}ms_start{startevt}end{endevt}_PoisToyData2D.csv")
