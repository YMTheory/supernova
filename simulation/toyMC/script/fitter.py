import numpy as np
import multiprocessing
from multiprocessing import cpu_count
from channel_analyser import channel


def scanning_asimov1D_separate(dt_arr, channels, MO):
    """
    1D NLL scanning: no C14 -> if C14, should use 2D fitting by default.
    :param dt_arr:
    :param channels:
    :param MO:
    :return:
    """
    nll = np.zeros_like(dt_arr)
    for idx, dT in enumerate(dt_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                val += cha.calc_Asimov_NLL_IO(dT, "MO")
            elif MO == "IO":
                val += cha.calc_Asimov_NLL_NO(dT, "MO")
        nll[idx] = val
    dchi2 = 2 * nll
    return dchi2


def scanning_asimov2D_separate(dt_arr, channels, MO):
    """
    2D NLL scanning, c14 can be included. Multi-processing is used here.
    :param dt_arr:
    :param channels:
    :param MO:
    :return:
    """
    nll = np.zeros_like(dt_arr)
    for idx, dt in enumerate(dt_arr):
        val = 0
        for cha in channels:
            if cha.name == "eES" or cha.name == "IBD":
                if MO == "NO":
                    val += cha.calc_Asimov_NLL_NO(dt, "MO")
                else:
                    val += cha.calc_Asimov_NLL_IO(dt, "MO")
        nll[idx] = val

    for cha in channels:
        if cha.name == "pES":
            if MO == "NO":
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nllpES = pool.map(cha.calc_Asimov_NLL_NO2D_withBkg, dt_arr)
            else:
                with multiprocessing.Pool(processes=cpu_count()) as pool:
                    nllpES = pool.map(cha.calc_Asimov_NLL_IO2D_withBkg, dt_arr)

        nll = nll + nllpES

    return nll


def scanning_asimov1D_combined(dt_arr, channels, MO, scaleIBD=1.0, tmin=-0.02, tmax=0.02, stept=0.00001):
    """
    combining different channels, as IBD and eES can not be separated in HyperK...
    :param dt_arr:
    :param cha_eES, cha_IBD
    :param MO:
    :return:
    """
    cha_eES, cha_IBD = None, None
    for cha in channels:
        if cha.name == "eES":
            cha_eES = cha
        elif cha.name == "IBD":
            cha_IBD = cha
    obs_t = np.arange(tmin, tmax, stept)
    if MO == "NO":
        obs_n = cha_eES._pdfNO_func(obs_t) * cha_eES.scale * stept + cha_IBD._pdfNO_func(obs_t) * scaleIBD * stept * cha_IBD.scale
    else:
        obs_n = cha_eES._pdfIO_func(obs_t) * cha_eES.scale * stept + cha_IBD._pdfIO_func(obs_t) * scaleIBD * stept * cha_IBD.scale

    nll = np.zeros_like(dt_arr)
    for idx, dt in enumerate(dt_arr):
        if MO == "NO":
            exp_n = cha_eES._pdfIO_func(obs_t + dt) * cha_eES.scale * stept + cha_IBD._pdfIO_func(obs_t + dt) * scaleIBD * stept * cha_IBD.scale
        else:
            exp_n = cha_eES._pdfNO_func(obs_t + dt) * cha_eES.scale * stept + cha_IBD._pdfNO_func(obs_t + dt) * scaleIBD * stept * cha_IBD.scale

        tmpnll = np.sum(obs_n - exp_n + obs_n * np.log(exp_n / obs_n))
        nll[idx] = tmpnll
        dchi2 = -2 * nll

    return dchi2

def scanning_asimov1D_partlysep(dt_arr, channels, MO, tagIBD=1.0, tmin=-0.02, tmax=0.02, stept=0.00001):
    dchi2_combined = scanning_asimov1D_combined(dt_arr, channels, MO, scaleIBD=1-tagIBD)
    for cha in channels:
        if cha.name == "IBD":
            dchi2_partIBD = scanning_asimov1D_separate(dt_arr, [cha], MO) * tagIBD
    dchi2 = dchi2_combined + dchi2_partIBD
    return dchi2


def find_locMin(dt_arr, dchi2_arr):
    idx = np.argmin(dchi2_arr)
    Tbest = dt_arr[idx]
    locMin = dchi2_arr[idx]
    return Tbest, locMin, idx

def generate_fine_dtarr(Tbest):
    TMIN = -0.009
    if Tbest < TMIN:
        Tbest = TMIN
    HALF_FITTING_NUM = 10
    TSTEP = 0.0001
    dt_arr = np.zeros(2*HALF_FITTING_NUM + 1)
    for i in range(2*HALF_FITTING_NUM + 1):
        dt_arr[i] = Tbest - HALF_FITTING_NUM * TSTEP + i * TSTEP
    return dt_arr

def parabola_fit(dt_arr, dchi2_arr):
    Tbest, locMin, idx = find_locMin(dt_arr, dchi2_arr)
    N = len(dt_arr)
    if idx < 2 or idx > N-3:
        print("**** Warning: best fit at the window edge!!! ****")
        return Tbest, locMin
    else:
        a, b, c = np.polyfit(dt_arr[idx-2:idx+3], dchi2_arr[idx-2:idx+3], 2)
        Tbest = - b / 2 /a
        locMin = (4*a*c - b**2) / 4 / a
        return Tbest, locMin, a, b, c


def scanning_chain(mode, dt_arr, channels, MO, tag_eff=0.5, tmin=-0.02, tmax=0.02, stept=0.00001):
    if mode == "1D_sep":
        dchi2_arr = scanning_asimov1D_separate(dt_arr, channels, MO)
    elif mode == "1D_comb":
        dchi2_arr = scanning_asimov1D_combined(dt_arr, channels, MO, tmin=tmin, tmax=tmax, stept=stept)
    elif mode == "1D_part":
        dchi2_arr = scanning_asimov1D_partlysep(dt_arr, channels, MO, tagIBD=tag_eff, tmin=tmin, tmax=tmax, stept=stept)
    Tbest, _, _ = find_locMin(dt_arr, dchi2_arr)
    dt_arr = generate_fine_dtarr(Tbest)
    if mode == "1D_sep":
        dchi2_arr = scanning_asimov1D_separate(dt_arr, channels, MO)
    elif mode == "1D_comb":
        dchi2_arr = scanning_asimov1D_combined(dt_arr, channels, MO, tmin=tmin, tmax=tmax, stept=stept)
    elif mode == "1D_part":
        dchi2_arr = scanning_asimov1D_partlysep(dt_arr, channels, MO, tagIBD=tag_eff, tmin=tmin, tmax=tmax, stept=stept)
    Tbest, locMin, a, b, c = parabola_fit(dt_arr, dchi2_arr)

    return dt_arr, dchi2_arr, Tbest, locMin, a, b, c

