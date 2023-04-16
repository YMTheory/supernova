import numpy as np
import ROOT

import warnings
warnings.filterwarnings("ignore")

import multiprocessing

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
            if cha == "pES":   # then do 2D fit
                dataT, dataE = cha.get_one_event2D(ievt)
                if MO == "NO":
                    val += cha.calc_NLL_NO2D(dataT, dataE, dT)
                else:
                    val += cha.calc_NLL_IO2D(dataT, dataE, dT)
            else:
                dataT = cha.get_one_event(ievt)
                if MO == "NO":
                    val += cha.calc_NLL_NO(dataT, dT)
                else:
                    val += cha.calc_NLL_IO(dataT, dT)
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
                    tmpval = cha.calc_Asimov_NLL_NO2D_withBkg(dT)
                    val += tmpval
                    #print("NO", cha.name, dT, tmpval)
                else:
                    tmpval = cha.calc_Asimov_NLL_NO(dT, "MO")
                    #print("NO", cha.name, dT, tmpval)
                    val += tmpval
            else:
                if cha.name == "pES":
                    tmpval = cha.calc_Asimov_NLL_IO2D_withBkg(dT)
                    #print("IO", cha.name, dT, tmpval)
                    val += tmpval
                else:
                    tmpval = cha.calc_Asimov_NLL_IO(dT, "MO")
                    #print("IO", cha.name, dT, tmpval)
                    val += tmpval

        #print(idx, dT, val)
        nll[idx] = val
    return nll


def scanning_asimov2D_withBkg_MP(dT_arr, channels, MO, ty):
    start_time = time.time()
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            #if cha.name == "IBD":
            if cha.name == "IBD" or cha.name == "eES":
                if MO == "NO":
                    tmpval = cha.calc_Asimov_NLL_NO(dT, "MO")
                else:
                    tmpval = cha.calc_Asimov_NLL_IO(dT, "MO")
                print(cha.name, dT, tmpval)
                val += tmpval
        nll[idx] = val

    for cha in channels:
        if cha.name == "pES":
            if MO == "NO":
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nllpES = pool.map(cha.calc_Asimov_NLL_NO2D_withBkg, dT_arr)
            else:
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nllpES = pool.map(cha.calc_Asimov_NLL_IO2D_withBkg, dT_arr)
            print("pES", nllpES)

            for i in range(len(dT_arr)):
                nll[i] = nll[i] + nllpES[i]
    

    #for cha in channels:
    #    if cha.name == "eES":
    #        if MO == "NO":
    #            with multiprocessing.Pool(processes=cpu_count()) as pool:
    #                nlleES = pool.map(cha.calc_Asimov_NLL_NO2D, dT_arr)
    #        else:
    #            with multiprocessing.Pool(processes=cpu_count()) as pool:
    #                nlleES = pool.map(cha.calc_Asimov_NLL_IO2D, dT_arr)
    #        print("eES", nlleES)
    #        
    #        for i in range(len(dT_arr)):
    #            nll[i] = nll[i] + nlleES[i]

    stop_time = time.time()
    print(f"Time collapsed for scanning_asimov2D_withBkg_MP is {stop_time - start_time}.")
    return nll


def scanning_asimov2D_MP(dT_arr, channels, MO, ty):
    start_time = time.time()
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            if cha.name == "IBD":
                if MO == "NO":
                    tmpval = cha.calc_Asimov_NLL_NO(dT, "MO")
                else:
                    tmpval = cha.calc_Asimov_NLL_IO(dT, "MO")

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

    for cha in channels:
        if cha.name == "eES":
            if MO == "NO":
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nlleES = pool.map(cha.calc_Asimov_NLL_NO2D, dT_arr)
            else:
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nlleES = pool.map(cha.calc_Asimov_NLL_IO2D, dT_arr)


    for i in range(len(dT_arr)):
        nll[i] = nll[i] + nlleES[i]
    stop_time = time.time()
    print(f"Time collapsed for scanning_asimov2D_MP is {stop_time - start_time}.")

    return nll


def scanning_asimov2D(dT_arr, channels, MO, ty):
    nll = np.zeros(len(dT_arr))
    for idx, dT in enumerate(dT_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                if cha.name == "pES":
                    tmpval = cha.calc_Asimov_NLL_NO2D(dT, ty=ty)
                    val += tmpval
                    #print(MO, cha.name, dT, tmpval)
                else:
                    tmpval = cha.calc_Asimov_NLL_NO(dT, ty)
                    val += tmpval
                    #print(MO, cha.name, dT, tmpval)

            if MO == "IO":
                if cha.name == "pES":
                    tmpval = cha.calc_Asimov_NLL_IO2D(dT, ty=ty)
                    val += tmpval
                    #print(MO, cha.name, dT, tmpval)
                else:
                    tmpval = cha.calc_Asimov_NLL_IO(dT, ty)
                    val += tmpval
                    #print(MO, cha.name, dT, tmpval)

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
    f = ROOT.TFile("/junofs/users/miaoyu/supernova/production/PDFs/backgrounds/C14/C14_rate_JUNO.root", "read")
    glow = f.Get("c14_low")
    return glow
    #ghig = f.Get("c14_high")
    #return glow, ghig






