import numpy as np
import multiprocessing
from multiprocessing import cpu_count
from channel_analyser import channel
from tqdm import tqdm
import uproot as up

def scanning_asimov1D_randomModel(dt_arr, channels, MO):
    nll = np.zeros_like(dt_arr)
    for i, dT in enumerate(dt_arr):
        val = 0
        for cha in channels:
            if MO == "NO":
                val += cha.calc_Asimov_NLL_NO(dT, "MO")
            else:
                val += cha.calc_Asimov_NLL_IO(dT, "MO")
        nll[i] = val
        dchi2 = 2 * nll

    return dchi2
def scanning_asimov1D_separate(dt_arr, channels, MO, reverse=False):
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
                if not reverse:
                    val += cha.calc_Asimov_NLL_IO(dT, "MO")
                else:
                    val += cha.calc_Asimov_NLL_NO(dT, "MO")
            elif MO == "IO":
                if not reverse:
                    val += cha.calc_Asimov_NLL_NO(dT, "MO")
                else:
                    val += cha.calc_Asimov_NLL_IO(dT, "MO")
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


def scanning_chain(mode, dt_arr, channels, MO, tag_eff=0.5, tmin=-0.02, tmax=0.02, stept=0.00001, reverse=False):
    if mode == "1D_sep":
        dchi2_arr = scanning_asimov1D_separate(dt_arr, channels, MO, reverse=reverse)
    elif mode == "1D_comb":
        dchi2_arr = scanning_asimov1D_combined(dt_arr, channels, MO, tmin=tmin, tmax=tmax, stept=stept)
    elif mode == "1D_part":
        dchi2_arr = scanning_asimov1D_partlysep(dt_arr, channels, MO, tagIBD=tag_eff, tmin=tmin, tmax=tmax, stept=stept)
    Tbest, _, _ = find_locMin(dt_arr, dchi2_arr)
    dt_arr = generate_fine_dtarr(Tbest)
    if mode == "1D_sep":
        dchi2_arr = scanning_asimov1D_separate(dt_arr, channels, MO, reverse=reverse)
    elif mode == "1D_comb":
        dchi2_arr = scanning_asimov1D_combined(dt_arr, channels, MO, tmin=tmin, tmax=tmax, stept=stept)
    elif mode == "1D_part":
        dchi2_arr = scanning_asimov1D_partlysep(dt_arr, channels, MO, tagIBD=tag_eff, tmin=tmin, tmax=tmax, stept=stept)
    Tbest, locMin, a, b, c = parabola_fit(dt_arr, dchi2_arr)

    return dt_arr, dchi2_arr, Tbest, locMin, a, b, c

def scanning_toyMC_chain(dt_arr, channels, MO, startid, endid):
    chi2 = np.zeros(endid-startid+1)
    for ievt in tqdm(range(startid, endid, 1)):
        tmpchi2 = scanning_toyMC_oneEvt(dt_arr, channels, ievt, MO)
        Tbest, _, _ = find_locMin(dt_arr, tmpchi2)
        fine_dt_arr = generate_fine_dtarr(Tbest)
        tmpchi2 = scanning_toyMC_oneEvt(fine_dt_arr, channels, ievt, MO)
        Tbest, locMin, a, b, c = parabola_fit(fine_dt_arr, tmpchi2)
        chi2[ievt-startid] = locMin
    return chi2

import random
def scanning_randomModels_chain(channels, MO, mode="random", dataModel="unknown", pdfModel="unknown", group="Garching"):
    Garching_models = [81123, 82503, 82703, 84003, 91123, 92503, 92703, 94003]
    if mode == "random":
        dataModel, pdfModel = random.sample(Garching_models, 2)
    print(dataModel, pdfModel)
    for cha in channels:
        cha.setMO(MO)
        cha.data_model = group
        cha.data_modelNo = dataModel
        cha.pdf_model = group
        cha.pdf_modelNo = pdfModel
        cha.setNOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{group}{pdfModel}_PDF_NO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
        cha.setIOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{group}{pdfModel}_PDF_IO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
        cha._load_pdf()
        cha.setNODataPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{group}{dataModel}_PDF_NO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
        cha.setIODataPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{group}{dataModel}_PDF_IO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
        cha._load_datapdf()

    dt_arrNO = np.arange(-0.010, 0.011, 0.001)
    chi2NO = scanning_asimov1D_randomModel(dt_arrNO, channels, "NO")
    Tbest, _, _, = find_locMin(dt_arrNO, chi2NO)
    dt_arrNO_fine = generate_fine_dtarr(Tbest)
    chi2NO = scanning_asimov1D_randomModel(dt_arrNO_fine, channels, "NO")
    TbestNO, locMinNO, aNO, bNO, cNO = parabola_fit(dt_arrNO_fine, chi2NO)

    dt_arrIO = np.arange(-0.010, 0.011, 0.001)
    chi2IO = scanning_asimov1D_randomModel(dt_arrIO, channels, "IO")
    Tbest, _, _, = find_locMin(dt_arrIO, chi2IO)
    dt_arrIO_fine = generate_fine_dtarr(Tbest)
    chi2IO = scanning_asimov1D_randomModel(dt_arrIO_fine, channels, "IO")
    TbestIO, locMinIO, aIO, bIO, cIO = parabola_fit(dt_arrIO_fine, chi2IO)

    print(locMinNO, locMinIO)
    if MO == "NO":
        dchi2 = locMinIO - locMinNO
    else:
        dchi2 = locMinNO - locMinIO
    return dchi2



def scanning_toyMC_oneEvt(dt_arr, channels, evtNO, pdfMO, reverse=False):
    nll = np.zeros(len(dt_arr))
    for i, dT in enumerate(dt_arr):
        for cha in channels:
            if cha.name == "eES" or cha.name == "IBD":
                oneevt = cha.get_one_event(evtNO)
                if pdfMO == "NO":
                    tmpnll = cha.calc_NLL_NO(oneevt, dT)
                else:
                    tmpnll = cha.calc_NLL_IO(oneevt, dT)
            elif cha.name == "pES":
                dataT, dataE = cha.get_one_event2D(evtNO)
                if pdfMO == "NO":
                    tmpnll = cha.calc_NLL_NO2D(dataT, dataE, dT)
                else:
                    tmpnll = cha.calc_NLL_IO2D(dataT, dataE, dT)

            nll[i] += tmpnll
    return 2 * nll


def gauss_fit(x, A, mu, sigma):
    return A * np.exp(-(x-mu)**2/2/sigma**2)

def load_C14(level):
    f = up.open("/junofs/users/miaoyu/supernova/production/PDFs/backgrounds/C14/C14_rate_JUNO.root")
    g = f["c14_"+level]

    x, y = g.values("x"), g.values("y")
    return x, y

def getC14Rate(level, Ethr):
    x, y = load_C14(level)
    return np.interp(Ethr, x, y)
