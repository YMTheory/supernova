import numpy as np
import ROOT
import time

import warnings
warnings.filterwarnings("ignore")

import multiprocessing
from multiprocessing import cpu_count

def scanning1D(dT_arr, channels, ievt, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
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
            if cha.name == "pES":   # then do 2D fit
                dataT, dataE = cha.get_one_event2D(ievt)
                if MO == "NO":
                    val += cha.calc_NLL_NO2D(dataT, dataE, dT)[0]
                else:
                    val += cha.calc_NLL_IO2D(dataT, dataE, dT)[0]
            else:
                dataT = cha.get_one_event(ievt)
                if MO == "NO":
                    val += cha.calc_NLL_NO(dataT, dT)
                else:
                    val += cha.calc_NLL_IO(dataT, dT)
            nll[idx] = val
    return nll


def scanning2D_absoluteMass(dT_arr, channels, ievt, MO):
    dt_used, nll = [], []
    for idx, dT in enumerate(dT_arr):
        val = 0
        totflag = False
        for cha in channels:
            dataT, dataE = cha.get_one_event2D(ievt)
            if MO == "NO":
                #val += cha.calc_NLL_NO2D(dataT, dataE, dT)
                tmp_val, flag = cha.calc_NLL_NO2D(dataT, dataE, dT)
                val += tmp_val
                if flag:
                    totflag = True
            else:
                #val += cha.calc_NLL_IO2D(dataT, dataE, dT)
                tmp_val, flag = cha.calc_NLL_IO2D(dataT, dataE, dT)
                val += tmp_val
                if flag:
                    totflag = True
        if not totflag:
            dt_used.append(dT)
            nll.append(val)

    dt_used = np.array(dt_used)
    nll = np.array(nll)
    return dt_used, nll


def scanning_asimov2D_allchannels(dT_arr, channels, ievt, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            dataT, dataE = cha.get_one_event2D(ievt)
            if MO == "NO":
                val += cha.calc_Asimov_NLL_NO2D(dataT, dataE, dT)
            else:
                val += cha.calc_Asimov_NLL_IO2D(dataT, dataE, dT)
        nll[idx] = val
    return nll


def scanning_asimov1D(dT_arr, channels, MO):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                tmpval = cha.calc_Asimov_NLL_NO(dT)
                #print(cha.name, dT, tmpval)
                val += tmpval
            else:
                tmpval = cha.calc_Asimov_NLL_IO(dT)
                #print(cha.name, dT, tmpval)
                val += tmpval
        nll[idx] = val

    return nll


def scanning_asimov2D_MP(dT_arr, channels, MO):
    start_time = time.time()
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            #if cha.name == "IBD":
            if cha.name == "IBD" or cha.name == "eES":
                if MO == "NO":
                    tmpval = cha.calc_Asimov_NLL_NO(dT)
                else:
                    tmpval = cha.calc_Asimov_NLL_IO(dT)
                val += tmpval
        nll[idx] = val

    for cha in channels:
        if cha.name == "pES":
            if MO == "NO":
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nllpES = pool.map(cha.calc_Asimov_NLL_NO2D, dT_arr)
            else:
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nllpES = pool.map(cha.calc_Asimov_NLL_IO2D, dT_arr)

            for i in range(len(dT_arr)):
                nll[i] = nll[i] + nllpES[i]
    

    stop_time = time.time()
    print(f"Time collapsed for scanning_asimov2D_MP is {stop_time - start_time}.")
    return nll


def scanning_asimov2D_allchannels_MP(dT_arr, channels, MO):
    start_time = time.time()
    nll = np.zeros(len(dT_arr))
    for cha in channels:
        #if cha.name == "eES" :
        if MO == "NO":
            with multiprocessing.Pool(processes=cpu_count()) as pool:
                tmpnll = pool.map(cha.calc_Asimov_NLL_NO2D, dT_arr)
        else:
            with multiprocessing.Pool(processes=cpu_count()) as pool:
                tmpnll = pool.map(cha.calc_Asimov_NLL_IO2D, dT_arr)
        for i in range(len(dT_arr)):
            nll[i] = nll[i] + tmpnll[i]

        #elif cha.name == "pES":
        #    nllpES = scanning_asimov1D(dT_arr, [cha], MO)
        #    for i in range(len(dT_arr)):
        #        nll[i] = nll[i] +  nllpES[i]

    stop_time = time.time()
    print(f"Time collapsed for scanning_asimov2D_allchannels_MP is {stop_time - start_time}.")

    return nll


def sanity_check(nll_arr):
    idx = np.argmin(nll_arr)
    if idx == len(nll_arr) - 1:
        return False
    else:
        return True
    
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
        print(f"Local minimum at {Tbest} with NLL = {locMin}.")
        print(f"Local minimum found at the edge at {dT_arr[idx]}.")
        if not param:
            return Tbest, locMin
        else:
            return Tbest, locMin, 0, 0, 0
    else:
        a, b, c = np.polyfit(dT_arr[idx-2:idx+3], nll_arr[idx-2:idx+3], 2)
        Tbest = - b / 2 / a
        locMin = (4*a*c - b**2) / 4 / a
        if not param:
            return Tbest, locMin
        else:
            return Tbest, locMin, a, b, c


def load_C14():
    f = ROOT.TFile("/junofs/users/miaoyu/supernova/production/PDFs/backgrounds/C14/C14_rate_JUNO.root", "read")
    ghig = f.Get("c14_high")
    glow = f.Get("c14_low")
    return glow, ghig



def scanning_asimov_chain(channels, MO, fitDim, param=False, details=False):
    ## NO Pdf fitting chain
    dt_arr = np.arange(-0.01, 0.011, 0.001)
    if fitDim == 2:
        nllNO_coarse = scanning_asimov2D_MP(dt_arr, channels, "NO")
    else:
        nllNO_coarse = scanning_asimov1D(dt_arr, channels, "NO")
    Tbest, locMin, _ = find_locMin(dt_arr, nllNO_coarse)
    dt_arr_fineNO = generate_fine_dTarr(Tbest)
    if fitDim == 2:
        nllNO_fine = scanning_asimov2D_MP(dt_arr_fineNO, channels, "NO")
    else:
        nllNO_fine = scanning_asimov1D(dt_arr_fineNO, channels, "NO" )
    TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dt_arr_fineNO, nllNO_fine, param=True)

    ## IO Pdf fitting chain
    dt_arr = np.arange(-0.01, 0.011, 0.001)
    if fitDim == 2:
        nllIO_coarse = scanning_asimov2D_MP(dt_arr, channels, "IO")
    else:
        nllIO_coarse = scanning_asimov1D(dt_arr, channels, "IO")
    Tbest, locMin, _ = find_locMin(dt_arr, nllIO_coarse)
    dt_arr_fineIO = generate_fine_dTarr(Tbest)
    if fitDim == 2:
        nllIO_fine = scanning_asimov2D_MP(dt_arr_fineIO, channels, "IO")
    else:
        nllIO_fine = scanning_asimov1D(dt_arr_fineIO, channels, "IO" )
    TbestFitIO, locMinFitIO, aIO, bIO, cIO = parabola_fit(dt_arr_fineIO, nllIO_fine, param=True)

    if MO == "NO":
        dchi2 = 2 * locMinFitIO - 2 * locMinFitNO
    else:
        dchi2 = 2 * locMinFitNO - 2 * locMinFitIO

    if not details:
        if not param:
            return TbestFitNO, locMinFitNO, TbestFitIO, locMinFitIO, dchi2
        else:
            return TbestFitNO, locMinFitNO, TbestFitIO, locMinFitIO, dchi2, aNO, bNO, cNO, aIO, bIO, cIO
    else:
        return dt_arr_fineNO, nllNO_fine, TbestFitNO, locMinFitNO, dt_arr_fineIO, nllIO_fine, TbestFitIO, locMinFitIO, dchi2, aNO, bNO, cNO, aIO, bIO, cIO


def scanning_asimov_chain_absoluteMass(channels, MO, fitDim):
    if MO == "NO":
        dt_arr = np.arange(-0.01, 0.011, 0.001)
        if fitDim == 2:
            nllNO_coarse = scanning_asimov2D_allchannels_MP(dt_arr, channels, MO)
        if fitDim == 1:
            nllNO_coarse = scanning_asimov1D(dt_arr, channels, MO)
        Tbest, locMin, _ = find_locMin(dt_arr, nllNO_coarse)
        dt_arr_fine = generate_fine_dTarr(Tbest)
        if fitDim == 2:
            nllNO_fine = scanning_asimov2D_allchannels_MP(dt_arr_fine, channels, MO)
        if fitDim == 1:
            nllNO_fine = scanning_asimov1D(dt_arr_fine, channels, MO)
        TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dt_arr_fine, nllNO_fine, param=True)
        return dt_arr_fine, nllNO_fine, TbestFitNO, locMinFitNO, aNO, bNO, cNO

    elif MO == "IO":
        dt_arr = np.arange(-0.01, 0.011, 0.001)
        if fitDim == 2:
            nllIO_coarse = scanning_asimov2D_allchannels_MP(dt_arr, channels, MO)
        if fitDim == 1:
            nllIO_coarse = scanning_asimov1D(dt_arr, channels, MO)
        Tbest, locMin, _ = find_locMin(dt_arr, nllIO_coarse)
        dt_arr_fine = generate_fine_dTarr(Tbest)
        if fitDim == 2:
            nllIO_fine = scanning_asimov2D_allchannels_MP(dt_arr_fine, channels, MO)
        if fitDim == 1:
            nllIO_fine = scanning_asimov1D(dt_arr_fine, channels, MO)
        TbestFitIO, locMinFitIO, aIO, bIO, cIO = parabola_fit(dt_arr_fine, nllIO_fine, param=True)

        return dt_arr_fine, nllIO_fine, TbestFitIO, locMinFitIO, aIO, bIO, cIO



def scanning_toyMC_chain_absoluteMass_2Dnew(channels, MO, evtNO):
    # directly scanning a large window -> avoid drop into local minimum.. Bad property of NLL
    dt_arr = np.arange(-0.01, 0.05, 0.001)
    if MO == "NO":
        dt_used, nll = scanning2D_absoluteMass(dt_arr, channels, evtNO, "NO")
        if len(dt_used) == 0:
            print("********** No time window satisties requirements.\n")
            return dt_used, nll, -10000, -10000, 0, 0, 0
        TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dt_used, nll, param=True)
        return dt_used, nll, TbestFitNO, locMinFitNO, aNO, bNO, cNO
    elif MO == "IO":
        dt_used, nll = scanning2D_absoluteMass(dt_arr, channels, evtNO, "IO")
        if len(dt_used) == 0:
            print("********** No time window satisties requirements.\n")
            return dt_used, nll, -10000, -10000, 0, 0, 0
        TbestFitIO, locMinFitIO, aIO, bIO, cIO = parabola_fit(dt_used, nll, param=True)
        return dt_used, nll, TbestFitIO, locMinFitIO, aIO, bIO, cIO



def scanning_toyMC_chain_absoluteMass(channels, MO, fitDim, evtNO):
    dt_arr = np.arange(-0.01, 0.011, 0.001)
    if MO == "NO":
        if fitDim == 2:
            nllNO_coarse = scanning2D_absoluteMass(dt_arr, channels, evtNO, "NO")
        Tbest, locMin, _ = find_locMin(dt_arr, nllNO_coarse)
        dt_arr_fine = generate_fine_dTarr(Tbest)
        if fitDim == 2:
            nllNO_fine = scanning2D_absoluteMass(dt_arr_fine, channels, evtNO, "NO")
        # sanity check :
        while True:
            if dt_arr_fine[-1] > 0.05:
                break
            if not sanity_check(nllNO_fine):
                print("Sanity check not satisfied -> run fine scanning again.")
                dt_arr_fine = generate_fine_dTarr(dt_arr_fine[-1])
                nllNO_fine = scanning2D_absoluteMass(dt_arr_fine, channels, evtNO, "NO")
            else:
                break
        TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dt_arr_fine, nllNO_fine, param=True)
        #print(f"Best fit point at {TbestFitNO} with {locMinFitNO}.") 
        return dt_arr_fine, nllNO_fine, TbestFitNO, locMinFitNO, aNO, bNO, cNO

    elif MO == "IO":
        if fitDim == 2:
            nllIO_coarse = scanning2D_absoluteMass(dt_arr, channels, evtNO, "IO")
        Tbest, locMin, _ = find_locMin(dt_arr, nllIO_coarse)
        dt_arr_fine = generate_fine_dTarr(Tbest)
        if fitDim == 2:
            nllIO_fine = scanning2D_absoluteMass(dt_arr_fine, channels, evtNO, "IO")
        # sanity check :
        while True:
            if dt_arr_fine[-1] > 0.05:
                break
            if not sanity_check(nllIO_fine):
                print("Sanity check not satisfied -> run fine scanning again.")
                dt_arr_fine = generate_fine_dTarr(dt_arr_fine[-1])
                nllIO_fine = scanning2D_absoluteMass(dt_arr_fine, channels, evtNO, "IO")
            else:
                break
            
        TbestFitIO, locMinFitIO, aIO, bIO, cIO = parabola_fit(dt_arr_fine, nllIO_fine, param=True)
        return dt_arr_fine, nllIO_fine, TbestFitIO, locMinFitIO, aIO, bIO, cIO



def scanning_toyMC_chain(channels, MO, fitDim, evtNO):
    ## NO Pdf fitting chain
    dt_arr = np.arange(-0.01, 0.011, 0.001)
    if fitDim == 2:
        nllNO_coarse = scanning2D(dt_arr, channels, evtNO,  "NO")
    else:
        nllNO_coarse = scanning1D(dt_arr, channels, evtNO, "NO")
    Tbest, locMin, _ = find_locMin(dt_arr, nllNO_coarse)
    dt_arr_fine = generate_fine_dTarr(Tbest)
    if fitDim == 2:
        nllNO_fine = scanning2D(dt_arr_fine, channels, evtNO, "NO")
    else:
        nllNO_fine = scanning1D(dt_arr, channels, evtNO, "NO" )
    TbestFitNO, locMinFitNO, aNO, bNO, cNO = parabola_fit(dt_arr_fine, nllNO_fine, param=True)

    ## IO Pdf fitting chain
    dt_arr = np.arange(-0.01, 0.011, 0.001)
    if fitDim == 2:
        nllIO_coarse = scanning2D(dt_arr, channels, evtNO, "IO")
    else:
        nllIO_coarse = scanning1D(dt_arr, channels, evtNO, "IO")
    Tbest, locMin, _ = find_locMin(dt_arr, nllIO_coarse)
    dt_arr_fine = generate_fine_dTarr(Tbest)
    if fitDim == 2:
        nllIO_fine = scanning2D(dt_arr_fine, channels, evtNO, "IO")
    else:
        nllIO_fine = scanning1D(dt_arr, channels, evtNO, "IO" )
    TbestFitIO, locMinFitIO, aIO, bIO, cIO = parabola_fit(dt_arr_fine, nllIO_fine, param=True)

    if MO == "NO":
        dchi2 = 2 * locMinFitIO - 2 * locMinFitNO
    else:
        dchi2 = 2 * locMinFitNO - 2 * locMinFitIO

    return TbestFitNO, locMinFitNO, TbestFitIO, locMinFitIO, dchi2



