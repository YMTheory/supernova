from curses.textpad import Textbox
import matplotlib as mpl
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

#plt.style.use("science")

import argparse

from tqdm import tqdm

from channel_analyser import channel


def plot_nll_one_event(dT, NO, IO):
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(dT, NO, "o-", ms=5, lw=2, label="NO")
    ax.plot(dT, IO, "o-", ms=5, lw=2, label="IO")
    ax.set_xlabel(r"$\Delta T$ [ms]", fontsize=15)
    ax.set_ylabel(r"NLL", fontsize=15)
    ax.legend(prop={"size":15})
    plt.show() 


def plot_nll_bundles(dT, NO, IO):
    fig, ax = plt.subplots(figsize=(6, 5))
    for oneNO, oneIO in zip(NO, IO):
        ax.plot(dT, oneNO, "-", lw=2, label="NO")
        ax.plot(dT, oneIO, "--", lw=2, label="IO")
    ax.set_xlabel(r"$\Delta T$ [ms]", fontsize=15)
    ax.set_ylabel(r"NLL", fontsize=15)
    plt.show() 

def plot_nll_fine(dT, arr):
    fig, ax = plt.subplots(figsize=(6, 5))
    for dTone, arrone in zip(dT, arr):
        ax.plot(dTone, arrone, "-", lw=2)
    ax.set_xlabel(r"$\Delta T$ [ms]", fontsize=15)
    ax.set_ylabel(r"NLL", fontsize=15)
    plt.savefig("fineNLL.pdf")
    #plt.show() 


if __name__ == "__main__":

    MO = "NO"
    model = "Garching82503"
    Ethr = 0.15
    
    parser = argparse.ArgumentParser(description='Arguments of SNNu analyser.')
    parser.add_argument('--model', type=str, default='Garching82503', help='Model name of SNNu.')
    parser.add_argument('--MO', type=str, default="NO", help="Mass ordering for the dataset.")
    parser.add_argument('--Ethr' , type=float, default=0.20, help="Detection threshold for pES channel, unit MeV.")
    args = parser.parse_args()
    
    model   = args.model
    Ethr    = args.Ethr
    MO      = args.MO
    
    channels = {}
    channels["pES"] = channel("pES", MO, model, Ethr)
    channels["IBD"] = channel("IBD", MO, model, 0.20)
    channels["eES"] = channel("eES", MO, model, 0.20)
    

    # Initialization, data loading...
    for cha in channels.values():
        cha._load_data()
        cha._load_pdf()

    #fig, ax = plt.subplots(figsize=(6, 4))

    # caluculating NLL and sensitivity:
    ### Coarse scanning:
    FITTING_EVENT_NUM =  500 # the sample number to run...
    dT_arr  = np.arange(-10, 11, 1) # coarse scannig, the step size is the bin width of PDFs
    nllNO   = np.zeros((FITTING_EVENT_NUM, len(dT_arr)))
    nllIO   = np.zeros((FITTING_EVENT_NUM, len(dT_arr)))
    #abnormalID = [140, 490]
    for ievt in tqdm(range(FITTING_EVENT_NUM)):
        #if ievt not in abnormalID:
        #    continue
        for idx, dT in enumerate(dT_arr):
            val_NO, val_IO = 0, 0
            for cha in channels.values():
                data = cha.get_one_event(ievt)
                val_NO += cha.calc_NLL_NO(data, dT)
                val_IO += cha.calc_NLL_IO(data, dT)
            
            nllNO[ievt, idx] = val_NO
            nllIO[ievt, idx] = val_IO

    TbestNO  = dT_arr[np.argmin(nllNO, axis=1)]
    locMinNO = np.min(nllNO, axis=1)
    TbestIO  = dT_arr[np.argmin(nllIO, axis=1)]
    locMinIO = np.min(nllIO, axis=1)

    ### Fine scanning:
    step_size, half_step_num = 0.1, 10
    dTNO, dTIO = [], []
    nllNO   = np.zeros((FITTING_EVENT_NUM, half_step_num * 2))
    nllIO   = np.zeros((FITTING_EVENT_NUM, half_step_num * 2))
    for ievt in tqdm(range(FITTING_EVENT_NUM)):
        dT_arr_NO = np.arange(TbestNO[ievt] - half_step_num * step_size, TbestNO[ievt] + half_step_num * step_size, step_size)
        dTNO.append(dT_arr_NO)
        dT_arr_IO = np.arange(TbestIO[ievt] - half_step_num * step_size, TbestIO[ievt] + half_step_num * step_size, step_size)
        dTIO.append(dT_arr_IO)
        for idx, dT in enumerate(dT_arr_NO):
            val_NO = 0
            for cha in channels.values():
                data = cha.get_one_event(ievt)
                val_NO += cha.calc_NLL_NO(data, dT)
            nllNO[ievt, idx] = val_NO
        for idx, dT in enumerate(dT_arr_IO):
            val_IO = 0
            for cha in channels.values():
                data = cha.get_one_event(ievt)
                val_IO += cha.calc_NLL_IO(data, dT)
            nllIO[ievt, idx] = val_IO

        #if ievt in [88, 255, 413]:
        #    ax.plot(dT_arr_IO, nllIO[ievt], ":", lw=2)

    dTNO = np.array(dTNO)
    dTIO = np.array(dTIO)

    ### parabola fitting :
    TbestNO, TbestIO = [], []
    locMinNO, locMinIO = [], []
    FITTING_EDGE_POINT = 2
    for ievt in range(FITTING_EVENT_NUM):
        #if ievt not in abnormalID:
        #    continue
        # check if the local minimum has enough neighbours for fitting.
        tmp_loc_min = np.argmin(nllNO[ievt])
        if tmp_loc_min < 2 or tmp_loc_min > half_step_num * 2 - 3:  ## no enough neighbours, not do fitting
            TbestNO.append(dTNO[ievt, tmp_loc_min])
            locMinNO.append(nllNO[ievt, tmp_loc_min])
        a, b, c = np.polyfit(dTNO[ievt, tmp_loc_min-2:tmp_loc_min+3], nllNO[ievt, tmp_loc_min-2:tmp_loc_min+3], 2)
        TbestNO.append(-b/2/a)
        locMinNO.append((4*a*c - b**2)/4/a)

        #### Abnormal logs
        if TbestNO[-1] <= -10 or TbestNO[-1] >= 10:
            print(f"One abnormal event {ievt}.")
            print(dTNO[ievt])
            print(nllNO[ievt])
            print(f"Fitting parameters {a}, {b}, {c}.")
            print(f"Fitted bestT = {TbestNO[-1]} and locMinNLL = {locMinNO[-1]}.")

        tmp_loc_min = np.argmin(nllIO[ievt])
        if tmp_loc_min < 2 or tmp_loc_min > half_step_num * 2 - 3:  ## no enough neighbours, not do fitting
            TbestIO.append(dTIO[ievt, tmp_loc_min])
            locMinIO.append(nllIO[ievt, tmp_loc_min])
        a, b, c = np.polyfit(dTIO[ievt, tmp_loc_min-2:tmp_loc_min+3], nllIO[ievt, tmp_loc_min-2:tmp_loc_min+3], 2)
        TbestIO.append(-b/2/a)
        locMinIO.append((4*a*c - b**2)/4/a)

        #### Abnormal logs
        if TbestIO[-1] <= -10 or TbestIO[-1] >= 10:
            print(f"One abnormal event {ievt}.")
            print(dTIO[ievt])
            print(nllIO[ievt])
            print(f"Fitting parameters {a}, {b}, {c}.")
            print(f"Fitted bestT = {TbestIO[-1]} and locMinNLL = {locMinIO[-1]}.")

        #if ievt in abnormalID:
        #    print(TbestIO[-1], locMinIO[-1])

    TbestNO = np.array(TbestNO)
    locMinNO = np.array(locMinNO)
    TbestIO = np.array(TbestIO)
    locMinIO = np.array(locMinIO)

    if MO == "NO" :
        sens = 2 * (locMinIO - locMinNO)
    else:
        sens = 2 * (locMinNO - locMinIO)


    #plt.savefig("check.pdf")
    #plt.show()

    ## plot single event NLL for checks
    # DISPLAY_EVENT_NUM = 10
    # plot_nll_one_event(dT_arr, nllNO[DISPLAY_EVENT_NUM], nllIO[DISPLAY_EVENT_NUM])
    # plot_nll_bundles(dT_arr, nllNO, nllIO)
    # plot_nll_fine(dTNO, nllNO)

    
    ### Save outputs from fitting
    df = pd.DataFrame({ 
                        "locMinNO" : locMinNO,
                        "locMinIO" : locMinIO,
                        "sens" : sens,
                        "TbestNO" : TbestNO,
                        "TbestIO" : TbestIO
                     })

    df.to_csv(f"./results/{model}_10kpc_{MO}_pESeESIBD_{Ethr:.2f}MeV.csv")
