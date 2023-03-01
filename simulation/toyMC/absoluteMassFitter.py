import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm

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


if __name__ == "__main__":
    MO          = "NO"
    model       = "Garching"
    modelNo     = 82703
    Ethr        = 0.20
    fitTmin     = -20
    fitTmax     = 30
    fileNo      = 0
    startevt    = 0
    endevt      = 0
    dataMass    = 0.0
    fitMass     = 0.0

    parser = argparse.ArgumentParser(description="ArgumentParser of the absolute mass fitter.")
    parser.add_argument('--model', type=str, default='Garching', help='Model name of SNNu.')
    parser.add_argument('--modelNo', type=int, default=82703, help="modelNo")
    parser.add_argument('--MO', type=str, default="NO", help="Mass ordering for the dataset.")
    parser.add_argument('--Ethr' , type=float, default=0.20, help="Detection threshold for pES channel, unit MeV.")
    parser.add_argument('--fitTmin', type=int, default=10, help="Minimum fitting time.")
    parser.add_argument('--fitTmax', type=int, default=50, help="Maximum fitting time.")
    parser.add_argument('--fileNo', type=int, default=0, help="data file number.")
    parser.add_argument("--start", type=int, default=0, help="Start event for fit.")
    parser.add_argument("--fitMass", type=float, default=0.0, help="Fit nuMass (eV).")
    parser.add_argument("--end", type=int, default=0, help="End event for fit.")
    args = parser.parse_args()

    MO          = args.MO
    model       = args.model
    modelNo     = args.modelNo
    fitTmin     = args.fitTmin
    fitTmax     = args.fitTmax
    fileNo      = args.fileNo
    startevt    = args.start
    endevt      = args.end
    fitMass     = args.fitMass

    channels = {}
    channels["eES"] = channel("eES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, nuMass=fitMass)   
    for cha in channels.values():
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
        nllNO_oneEvt = scanning2D(dT_arr, channels.values(), ievt-startevt, "NO")
        #print("Rough nllNO_oneEvt ", nllNO_oneEvt)
        Tbest, locMin, _ = find_locMin(dT_arr, nllNO_oneEvt)
        #print(Tbest, locMin)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        nllNO_oneEvt = scanning2D(dT_arr_fine, channels.values(), ievt-startevt, "NO")
        #print("Fine  nllNO_oneEvt ", nllNO_oneEvt)
        TbestFit, locMinFit = parabola_fit(dT_arr_fine, nllNO_oneEvt)
        #print(TbestFit, locMinFit)
    
        TbestNO[ievt-startevt] = TbestFit
        locMinNO[ievt-startevt] = locMinFit
    
    for ievt in tqdm(range(startevt, endevt, 1)):
        nllIO_oneEvt = scanning2D(dT_arr, channels.values(), ievt-startevt, "IO")
        Tbest, locMin, _ = find_locMin(dT_arr, nllIO_oneEvt)
        dT_arr_fine = generate_fine_dTarr(Tbest)
        nllIO_oneEvt = scanning2D(dT_arr_fine, channels.values(), ievt-startevt, "IO")
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

    df.to_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}{modelNo}_10kpc_{MO}_eESonly_{Ethr:.2f}MeV_{fitMass:.1f}eV_fitTmax{fitTmax:d}ms_start{startevt}end{endevt}_PoisToyData2D.csv")
    


