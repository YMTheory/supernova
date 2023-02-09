import sys
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy import special
from prettytable import PrettyTable

def read(model, cha, E, T, dist, end="tot"):

    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO = [], [], [], []
    sens_dataNO, sens_dataIO = [], []
    for i in range(1):
        #df1 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_NO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")
        #df2 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_IO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")
        df1 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_NO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")
        df2 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_IO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")

        tmp_deltaT_dataNO_pdfNO = df1["TbestNO"]
        tmp_deltaT_dataNO_pdfIO = df1["TbestIO"]
        tmp_sens_dataNO = df1["sens"]

        tmp_deltaT_dataIO_pdfNO = df2["TbestNO"]
        tmp_deltaT_dataIO_pdfIO = df2["TbestIO"]
        tmp_sens_dataIO = df2["sens"]

        for a, b, c, d, e, f in zip(tmp_deltaT_dataNO_pdfNO, tmp_deltaT_dataNO_pdfIO, tmp_sens_dataNO, tmp_deltaT_dataIO_pdfNO, tmp_deltaT_dataIO_pdfIO, tmp_sens_dataIO):
            deltaT_dataNO_pdfNO.append(a)
            deltaT_dataNO_pdfIO.append(b)
            sens_dataNO.append(c)
            deltaT_dataIO_pdfNO.append(d)
            deltaT_dataIO_pdfIO.append(e)
            sens_dataIO.append(f)

    deltaT_dataNO_pdfNO = np.array(deltaT_dataNO_pdfNO)
    deltaT_dataNO_pdfIO = np.array(deltaT_dataNO_pdfIO)
    sens_dataNO = np.array(sens_dataNO)
    deltaT_dataIO_pdfNO = np.array(deltaT_dataIO_pdfNO)
    deltaT_dataIO_pdfIO = np.array(deltaT_dataIO_pdfIO)
    sens_dataIO = np.array(sens_dataIO)

    return deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO 



def calcSens():
    model = "Garching82703"
    cha = "pESeESIBD"
    E = 0.10
    T = 20
    dist = 10
    _, _, sens_dataNO, _, _, sens_dataIO = read(model, cha, E, T, dist, end="tot_doFit_binnedData_unbinnedNLL_C14bkghigh_PosiToyData")

    medsens_NO = np.median(sens_dataNO)
    medsens_IO = np.median(sens_dataIO)

    nNOtest, nIOtest = 0, 0
    for i, j in zip(sens_dataNO, sens_dataIO):
        if i < -medsens_IO:
            nNOtest += 1
        if j < -medsens_NO:
            nIOtest += 1
    print(nNOtest, nIOtest)
    nNOtest /= len(sens_dataNO)
    nIOtest /= len(sens_dataIO)
    
    nsigma_NO = np.sqrt(2) * special.erfcinv(2*nNOtest)
    nsigma_IO = np.sqrt(2) * special.erfcinv(2*nIOtest)

    asimovNO_dchi2 = 7.3
    asimovIO_dchi2 = 8.4

    table = PrettyTable(['MO','dchi square','p-value','n sigma', 'Asimov sensitivity'])
    table.add_row(['NO', f"{medsens_NO:.3f}", nNOtest, f"{nsigma_NO:.2f}", f"{np.sqrt(asimovNO_dchi2):.2f}"])
    table.add_row(['IO', f"{medsens_IO:.3f}", nIOtest, f"{nsigma_IO:.2f}", f"{np.sqrt(asimovIO_dchi2):.2f}"])
    print(table)


if __name__ == "__main__" :
    calcSens()


