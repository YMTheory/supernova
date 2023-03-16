import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import pyplot as plt, cm
from matplotlib import colors
import matplotlib as mpl
plt.style.use("science")
import ROOT
import uproot as up
from datetime import datetime


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


def read1(filename1, filename2):

    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO = [], [], [], []
    sens_dataNO, sens_dataIO = [], []
    for i in range(1):
        #df1 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_NO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")
        #df2 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_IO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")
        df1 = pd.read_csv(filename1)
        df2 = pd.read_csv(filename2)

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

def read2(filename1, filename2):

    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO = [], [], [], []
    sens_dataNO, sens_dataIO = [], []
    NeES_NO, NpES_NO, NIBD_NO =  [], [], []
    NeES_IO, NpES_IO, NIBD_IO =  [], [], []
    for i in range(1):
        #df1 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_NO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")
        #df2 = pd.read_csv(f"/junofs/users/miaoyu/supernova/simulation/toyMC/results/{model}_{dist}kpc_IO_{cha}_{E:.2f}MeV_fitTmax{T}ms_{end}.csv")
        df1 = pd.read_csv(filename1)
        df2 = pd.read_csv(filename2)

        tmp_deltaT_dataNO_pdfNO = df1["TbestNO"]
        tmp_deltaT_dataNO_pdfIO = df1["TbestIO"]
        tmp_sens_dataNO = df1["sens"]
        tmp_NeES_NO = df1["NeES"]
        tmp_NpES_NO = df1["NpES"]
        tmp_NIBD_NO = df1["NIBD"]

        tmp_deltaT_dataIO_pdfNO = df2["TbestNO"]
        tmp_deltaT_dataIO_pdfIO = df2["TbestIO"]
        tmp_sens_dataIO = df2["sens"]
        tmp_NeES_IO = df2["NeES"]
        tmp_NpES_IO = df2["NpES"]
        tmp_NIBD_IO = df2["NIBD"]

        for a, b, c, d, e, f in zip(tmp_deltaT_dataNO_pdfNO, tmp_deltaT_dataNO_pdfIO, tmp_sens_dataNO, tmp_deltaT_dataIO_pdfNO, tmp_deltaT_dataIO_pdfIO, tmp_sens_dataIO):
            deltaT_dataNO_pdfNO.append(a)
            deltaT_dataNO_pdfIO.append(b)
            sens_dataNO.append(c)
            deltaT_dataIO_pdfNO.append(d)
            deltaT_dataIO_pdfIO.append(e)
            sens_dataIO.append(f)

        for n1, n2, n3, n4, n5, n6 in zip(tmp_NeES_NO, tmp_NpES_NO, tmp_NIBD_NO, tmp_NeES_IO, tmp_NpES_IO, tmp_NIBD_IO):
            NeES_NO.append(n1)
            NpES_NO.append(n2)
            NIBD_NO.append(n3)
            NeES_IO.append(n4)
            NpES_IO.append(n5)
            NIBD_IO.append(n6)

    deltaT_dataNO_pdfNO = np.array(deltaT_dataNO_pdfNO)
    deltaT_dataNO_pdfIO = np.array(deltaT_dataNO_pdfIO)
    sens_dataNO = np.array(sens_dataNO)
    deltaT_dataIO_pdfNO = np.array(deltaT_dataIO_pdfNO)
    deltaT_dataIO_pdfIO = np.array(deltaT_dataIO_pdfIO)
    sens_dataIO = np.array(sens_dataIO)
    NeES_NO = np.array(NeES_NO)
    NpES_NO = np.array(NpES_NO)
    NIBD_NO = np.array(NIBD_NO)
    NeES_IO = np.array(NeES_IO)
    NpES_IO = np.array(NpES_IO)
    NIBD_IO = np.array(NIBD_IO)

    return deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO, NeES_NO, NpES_NO, NIBD_NO, NeES_IO, NpES_IO, NIBD_IO

def read3(filename1, filename2):

    minNLL_dataNO_pdfNO, minNLL_dataNO_pdfIO, minNLL_dataIO_pdfNO, minNLL_dataIO_pdfIO = [], [], [], []
    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO = [], [], [], []
    for i in range(1):
        df1 = pd.read_csv(filename1)
        df2 = pd.read_csv(filename2)

        tmp_locMinNLL_dataNO_pdfNO = df1["locMinNO"]
        tmp_locMinNLL_dataNO_pdfIO = df1["locMinIO"]
        tmp_deltaT_dataNO_pdfNO = df1["TbestNO"]
        tmp_deltaT_dataNO_pdfIO = df1["TbestIO"]

        tmp_locMinNLL_dataIO_pdfNO = df2["locMinNO"]
        tmp_locMinNLL_dataIO_pdfIO = df2["locMinIO"]
        tmp_deltaT_dataIO_pdfNO = df1["TbestNO"]
        tmp_deltaT_dataIO_pdfIO = df1["TbestIO"]

        for a, b, c, d in zip(tmp_locMinNLL_dataNO_pdfNO, tmp_locMinNLL_dataNO_pdfIO, tmp_locMinNLL_dataIO_pdfNO, tmp_locMinNLL_dataIO_pdfIO):
            minNLL_dataNO_pdfNO.append(a)
            minNLL_dataNO_pdfIO.append(b)
            minNLL_dataIO_pdfNO.append(c)
            minNLL_dataIO_pdfIO.append(d)
        for a, b, c, d in zip(tmp_deltaT_dataNO_pdfNO, tmp_deltaT_dataNO_pdfIO, tmp_deltaT_dataIO_pdfNO, tmp_deltaT_dataIO_pdfIO):
            deltaT_dataNO_pdfNO.append(a)
            deltaT_dataNO_pdfIO.append(b)
            deltaT_dataIO_pdfNO.append(c)
            deltaT_dataIO_pdfIO.append(d)

    minNLL_dataNO_pdfNO = np.array(minNLL_dataNO_pdfNO)    
    minNLL_dataNO_pdfIO = np.array(minNLL_dataNO_pdfIO)    
    minNLL_dataIO_pdfNO = np.array(minNLL_dataIO_pdfNO)    
    minNLL_dataIO_pdfIO = np.array(minNLL_dataIO_pdfIO)    
    deltaT_dataNO_pdfNO = np.array(deltaT_dataNO_pdfNO)    
    deltaT_dataNO_pdfIO = np.array(deltaT_dataNO_pdfIO)    
    deltaT_dataIO_pdfNO = np.array(deltaT_dataIO_pdfNO)    
    deltaT_dataIO_pdfIO = np.array(deltaT_dataIO_pdfIO)    

    return minNLL_dataNO_pdfNO, minNLL_dataNO_pdfIO, minNLL_dataIO_pdfNO, minNLL_dataIO_pdfIO, deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO

def drawHist2d():
    minNLL_dataNO_pdfNO, minNLL_dataNO_pdfIO, minNLL_dataIO_pdfNO, minNLL_dataIO_pdfIO, deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO = read3("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_0311.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv")
    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO, NeES_NO, NpES_NO, NIBD_NO, NeES_IO, NpES_IO, NIBD_IO = read2("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_0311.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv")

    hNO = ROOT.TH2D("hNO", "", 100, 0, 0.01,  100, -30, 30)
    hIO = ROOT.TH2D("hIO", "", 100, -0.01, 0, 100, -30, 30)
    
    for i, j in zip(deltaT_dataNO_pdfNO-deltaT_dataNO_pdfIO, minNLL_dataNO_pdfNO-minNLL_dataNO_pdfIO):
        hNO.Fill(i, 2*j)

    for i, j in zip(deltaT_dataIO_pdfIO-deltaT_dataIO_pdfNO, minNLL_dataIO_pdfIO-minNLL_dataIO_pdfNO):
        hIO.Fill(i, 2*j)

    arrNO, arrIO = np.zeros((hNO.GetNbinsY(), hNO.GetNbinsX())), np.zeros((hIO.GetNbinsY(), hIO.GetNbinsX()))
    for i in range(hNO.GetNbinsY()):
        for j in range(hNO.GetNbinsX()):
            arrNO[i, j] = hNO.GetBinContent(j+1, hNO.GetNbinsY()-i-1)
    for i in range(hIO.GetNbinsY()):
        for j in range(hIO.GetNbinsX()):
            arrIO[i, j] = hIO.GetBinContent(j+1, hIO.GetNbinsY()-i-1)

    fig ,ax = plt.subplots(1, 2, figsize=(10, 4))

    im0 = ax[0].imshow(arrNO, extent=[0, 0.01, -30, 30], aspect="auto", norm=colors.LogNorm(vmin=1))
    cb0 = plt.colorbar(im0, ax=ax[0])                  
    ax[0].set_xlabel(r"$\Delta t$(NO) - $\Delta t$(IO) [s]", fontsize=18)   
    im1 = ax[1].imshow(arrIO, extent=[-0.01, 0, -30, 30], aspect="auto", norm=colors.LogNorm(vmin=1))
    cb1 = plt.colorbar(im1, ax=ax[1])
    ax[1].set_xlabel(r"$\Delta t$(IO) - $\Delta t$(NO) [s]", fontsize=18)   

    ax[0].set_ylabel(r"$\Delta \chi^2$ (NO-IO)", fontsize=18)
    ax[0].set_title("NO data", fontsize=18)
    ax[0].tick_params(axis="x", labelsize=18)
    ax[0].tick_params(axis="y", labelsize=18)

    ax[1].set_ylabel(r"$\Delta \chi^2$ (IO-NO)", fontsize=18)
    ax[1].set_title("IO data", fontsize=18)
    ax[1].tick_params(axis="x", labelsize=18)
    ax[1].tick_params(axis="y", labelsize=18)
    plt.tight_layout()
    currentDateAndTime = datetime.now()
    plt.savefig(f"../plots/Garching82703_pESeESIBD_10kpc_0.10MeV_noC14bkg_fit1DTobs_Tmin-20msTmax20ms_NLLdeltaTHist2D_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}.pdf")
    plt.show()


def scatter():
    minNLL_dataNO_pdfNO, minNLL_dataNO_pdfIO, minNLL_dataIO_pdfNO, minNLL_dataIO_pdfIO, deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO = read3("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv")
    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO, NeES_NO, NpES_NO, NIBD_NO, NeES_IO, NpES_IO, NIBD_IO = read2("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv")

    #fig ,ax = plt.subplots(1, 2, figsize=(10, 4))
    #ax[0].plot(deltaT_dataNO_pdfNO-deltaT_dataNO_pdfIO, minNLL_dataNO_pdfNO-minNLL_dataNO_pdfIO, "o", ms=1, alpha=0.5)
    #ax[0].set_xlabel(r"$\Delta t$ [s]", fontsize=18)
    #ax[1].plot(deltaT_dataIO_pdfNO-deltaT_dataIO_pdfIO, minNLL_dataIO_pdfNO-minNLL_dataIO_pdfIO, "o", ms=1, alpha=0.5)
    #ax[1].set_xlabel(r"$\Delta t$ [s]", fontsize=18)

    #ax[0].set_ylabel(r"$\Delta$ NLL", fontsize=18)
    #ax[0].set_title("NO data", fontsize=18)
    #ax[0].tick_params(axis="x", labelsize=18)
    #ax[0].tick_params(axis="y", labelsize=18)

    #ax[1].set_ylabel(r"$\Delta$ NLL", fontsize=18)
    #ax[1].set_title("IO data", fontsize=18)
    #ax[1].tick_params(axis="x", labelsize=18)
    #ax[1].tick_params(axis="y", labelsize=18)


    fig ,ax = plt.subplots(figsize=(7, 4))
    ax.plot(NeES_NO, minNLL_dataNO_pdfNO-minNLL_dataNO_pdfIO, "o", ms=2, alpha=0.5, fillstyle="none", label="NO data")
    ax.plot(NeES_IO, minNLL_dataIO_pdfNO-minNLL_dataIO_pdfIO, "o", ms=2, alpha=0.5, fillstyle="none", label="IO data")
    ax.set_xlabel("eES count", fontsize=18)
    ax.set_ylabel(r"$\Delta$ NLL", fontsize=18)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    ax.legend(prop={"size":14})


    plt.tight_layout()
    currentDateAndTime = datetime.now()
    plt.savefig(f"../plots/Garching82703_eESeESeES_10kpc_0.10MeV_noC14bkg_fit1DTobs_Tmin-20msTmax20ms_NLLNeESscatter_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}.pdf")
    plt.show()


def drawNLL():
    drawHist = False
    minNLL_dataNO_pdfNO, minNLL_dataNO_pdfIO, minNLL_dataIO_pdfNO, minNLL_dataIO_pdfIO, _, _, _, _ = read3("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv")

    fig ,ax = plt.subplots(1, 2, figsize=(10, 4))
    if drawHist:
        ax[0].hist(minNLL_dataNO_pdfNO, bins=100, range=(-500, 0), alpha=0.3, color="blue", label=f"NO PDF: {np.median(minNLL_dataNO_pdfNO):.4e}")
        ax[0].hist(minNLL_dataNO_pdfIO, bins=100, range=(-500, 0), alpha=0.3, color="red",  label=f"IO PDF: {np.median(minNLL_dataNO_pdfIO):.4e}")
        ax[0].set_xlabel(r"$\Delta t$ [ms]", fontsize=18)
        ax[0].set_title("NO data", fontsize=18)
        ax[0].tick_params(axis="x", labelsize=18)
        ax[0].tick_params(axis="y", labelsize=0)
        ax[0].legend(prop={"size":15})
        #ax[0].semilogy()

        ax[1].hist(minNLL_dataIO_pdfNO, bins=100, range=(-500, 0), alpha=0.3, color="blue", label=f"NO PDF: {np.median(minNLL_dataIO_pdfNO):.4e}")
        ax[1].hist(minNLL_dataIO_pdfIO, bins=100, range=(-500, 0), alpha=0.3, color="red",  label=f"IO PDF: {np.median(minNLL_dataIO_pdfIO):.4e}")
        ax[1].set_xlabel(r"$\Delta t$ [ms]", fontsize=18)
        ax[1].set_title("IO data", fontsize=18)
        ax[1].tick_params(axis="x", labelsize=18)
        ax[1].tick_params(axis="y", labelsize=0)
        ax[1].legend(prop={"size":15})
        #ax[1].semilogy()

    else:
        ax[0].plot(minNLL_dataNO_pdfNO, minNLL_dataNO_pdfIO, "o", ms=1, alpha=0.5)
        ax[0].plot([np.min(minNLL_dataNO_pdfNO), np.max(minNLL_dataNO_pdfNO)], [np.min(minNLL_dataNO_pdfNO), np.max(minNLL_dataNO_pdfNO)], "--", lw=2, color="red")
        ax[0].set_xlabel("NLL (NO fit)", fontsize=18)
        ax[0].set_ylabel("NLL (IO fit)", fontsize=18)
        ax[0].set_title("NO data", fontsize=18)
        ax[0].tick_params(axis="x", labelsize=18)
        ax[0].tick_params(axis="y", labelsize=18)

        ax[1].plot(minNLL_dataIO_pdfIO, minNLL_dataIO_pdfNO, "o", ms=1, alpha=0.5)
        ax[1].plot([np.min(minNLL_dataIO_pdfIO), np.max(minNLL_dataIO_pdfIO)], [np.min(minNLL_dataIO_pdfIO), np.max(minNLL_dataIO_pdfIO)], "--", lw=2, color="red")
        ax[1].set_xlabel("NLL (NO fit)", fontsize=18)
        ax[1].set_ylabel("NLL (IO fit)", fontsize=18)
        ax[1].set_title("IO data", fontsize=18)
        ax[1].tick_params(axis="x", labelsize=18)
        ax[1].tick_params(axis="y", labelsize=18)

    plt.tight_layout()
    currentDateAndTime = datetime.now()
    plt.savefig(f"../plots/Garching82703_eESIBD_10kpc_0.10MeV_noC14bkg_fit1DTobs_Tmin-20msTmax20ms_minLocNLLScat_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}.pdf")
    plt.show()




def drawFitT():

    #deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO = read1("../results/Garching82703_10kpc_NO_eESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D.csv", "../results/Garching82703_10kpc_IO_eESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D.csv")
    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO = read1("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_v2.csv")

    fig ,ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(deltaT_dataNO_pdfNO, bins=100, range=(-0.01, 0.01), alpha=0.3, color="blue", label=f"NO PDF: {np.mean(deltaT_dataNO_pdfNO)*1000:.2f} ms")
    ax[0].hist(deltaT_dataNO_pdfIO, bins=100, range=(-0.01, 0.01), alpha=0.3, color="red",  label=f"IO PDF: {np.mean(deltaT_dataNO_pdfIO)*1000:.2f} ms")
    ax[0].set_xlabel(r"$\Delta t$ [ms]", fontsize=18)
    ax[0].set_title("NO data", fontsize=18)
    ax[0].tick_params(axis="x", labelsize=18)
    ax[0].tick_params(axis="y", labelsize=0)
    ax[0].legend(prop={"size":15})
    ax[0].vlines(0, 0, 10000, linestyle="--", lw=2, color="black")
    #ax[0].semilogy()

    ax[1].hist(deltaT_dataIO_pdfNO, bins=100, range=(-0.01, 0.01), alpha=0.3, color="blue", label=f"NO PDF: {np.mean(deltaT_dataIO_pdfNO)*1000:.2f} ms")
    ax[1].hist(deltaT_dataIO_pdfIO, bins=100, range=(-0.01, 0.01), alpha=0.3, color="red",  label=f"IO PDF: {np.mean(deltaT_dataIO_pdfIO)*1000:.2f} ms")
    ax[1].set_xlabel(r"$\Delta t$ [ms]", fontsize=18)
    ax[1].set_title("IO data", fontsize=18)
    ax[1].tick_params(axis="x", labelsize=18)
    ax[1].tick_params(axis="y", labelsize=0)
    ax[1].legend(prop={"size":15})
    ax[1].vlines(0, 0, 10000, linestyle="--", lw=2, color="black")
    #ax[1].semilogy()

    plt.tight_layout()
    currentDateAndTime = datetime.now()
    plt.savefig(f"../plots/Garching82703_eESIBD_10kpc_0.10MeV_noC14bkg_fit1DTobs_Tmin-20msTmax20ms_deltaFitT_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}_v2.pdf")
    plt.show()


def draw():
    # Ethr = 0.15 MeV
    asimovNO = 7.20
    asimovIO = 7.59
    # Ethr = 0.10 MeV
    asimovNO = 20.3
    asimovIO = 14.7

    Et = 0.10
    fig, ax = plt.subplots(figsize=(9, 4))
    #_, _, sens_dataNO, _, _, sens_dataIO                = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyData2D")
    #deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO = read1("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D.csv")
    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO = read1("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.60MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_THEIA100.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.60MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_THEIA100.csv")
    h1 = ax.hist(sens_dataNO,   bins=100, range=(-50, 50),   alpha=0.5, label="true NO", color="blue",)
    h2 = ax.hist(-sens_dataIO,  bins=100, range=(-50, 50),  alpha=0.5, label="true IO", color="red")


    conts1, edges1 = np.histogram(sens_dataNO, bins=100, range=(-50, 50))
    conts2, edges2 = np.histogram(-sens_dataIO, bins=100, range=(-50, 50))
    cents1, cents2 = [], [],
    for i in range(100):
        cents1.append((edges1[i] + edges1[i+1])/2.)
        cents2.append((edges2[i] + edges2[i+1])/2.)

    def gauss(x, A, mu, sigma):
        return A * np.exp(-(x-mu)**2/2/sigma**2)

    from scipy.optimize import curve_fit
    popt1, pcov1 = curve_fit(gauss, cents1, conts1)
    print(popt1)
    popt2, pcov2 = curve_fit(gauss, cents2, conts2)
    print(popt2)
    dx1 = np.arange(edges1[0], edges1[-1], 0.1)
    dy1 = gauss(dx1, *popt1)
    dx2 = np.arange(edges2[0], edges2[-1], 0.1)
    dy2 = gauss(dx2, *popt2)
    ax.plot(dx1, dy1, "-", color="black", lw=2)
    ax.plot(dx2, dy2, "-", color="black", lw=2)
    ax.text(33, 2000, "gauss", fontsize=15)


    currentDateAndTime = datetime.now()
    
    ax.set_xlabel(r"$\Delta\chi^2$", fontsize=18)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=0)
    #ax.set_ylim(0, 7000)

    l1 = ax.vlines(-np.median(sens_dataIO), 0, 3500, linestyle=":", lw=2.5, color="darkviolet")
    l3 = ax.vlines(-asimovIO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(-np.median(sens_dataIO)-6, 1000, f"{-np.median(sens_dataIO):.1f}", fontsize=18, color="red")
    ax.text(-asimovIO+2, 1500, f"{-asimovIO:.1f}", fontsize=18, color="black")

    ax.vlines(np.median(sens_dataNO),  0, 3500, linestyle=":", lw=2.5, color="darkviolet")
    ax.vlines(asimovNO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(np.median(sens_dataNO)-6, 1000, f"{np.median(sens_dataNO):.1f}", fontsize=18, color="blue")
    ax.text(asimovNO+2, 1500, f"{asimovNO:.1f}", fontsize=18, color="black")

    
    ax.set_ylim(0, 7000)
    lg1 = ax.legend([l1, l3], ["median value of Monte Carlo", "Asimov dataset"], loc="upper right", prop={"size":16},frameon=True)
    plt.legend(prop={"size":16}, loc="upper left", frameon=True, ncol=3)
    plt.gca().add_artist(lg1)

    plt.tight_layout()
    plt.savefig(f"../plots/Garching82703_pESeESIBD_10kpc_0.60MeV_noC14bkg_fit1DTobs_Tmin-20msTmax20ms_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}_THEIA100.pdf")
    plt.show()



def drawC14():
    Et = 0.10
    _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    _, _, sens_dataNO_low, _, _, sens_dataIO_low    = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_C14bkglow_PosiToyData")
    _, _, sens_dataNO_high, _, _, sens_dataIO_high    = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_C14bkghigh_PosiToyData")

    fig, ax = plt.subplots(figsize=(9, 4))
    h1 = ax.hist(sens_dataNO,   bins=100, range=(-30, 30),   alpha=0.5, label="true NO (no C14)", color="blue",)
    h2 = ax.hist(-sens_dataIO,  bins=100, range=(-30, 30),  alpha=0.5, label="true IO (no C14)", color="red")
    h3 = ax.hist(sens_dataNO_low,   bins=100, range=(-30, 30),   alpha=0.5, label="true NO (1e-18 g/g)", histtype="step", lw=2, color="green",)
    h4 = ax.hist(-sens_dataIO_low,  bins=100, range=(-30, 30),  alpha=0.5, label="true IO (1e-18 g/g)", histtype="step", lw=2, color="orange")
    h5 = ax.hist(sens_dataNO_high,   bins=100, range=(-30, 30),  linestyle=":",  alpha=0.5, label="true NO (1e-17 g/g)", histtype="step", lw=2, color="darkviolet",)
    h6 = ax.hist(-sens_dataIO_high,  bins=100, range=(-30, 30),  linestyle=":", alpha=0.5, label="true IO (1e-17 g/g)", histtype="step", lw=2, color="brown")
    #h1 = ax.hist(sens_dataNO,   bins=100, range=(-30, 30),   alpha=0.5, label="true NO (unbinned)", color="blue",)
    #h2 = ax.hist(-sens_dataIO,  bins=100, range=(-30, 30),  alpha=0.5, label="true IO (unbinned)", color="red")
    #h3 = ax.hist(sens_dataNO_low,   bins=100, range=(-30, 30),   alpha=0.5, label="true NO (0.01 ms binning)", histtype="step", lw=2, color="green",)
    #h4 = ax.hist(-sens_dataIO_low,  bins=100, range=(-30, 30),  alpha=0.5, label="true IO (0.01 ms binning)", histtype="step", lw=2, color="orange")
    ax.set_xlabel(r"$\Delta\chi^2$", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=0)
    #ax.set_ylim(0, 7000)

    l1 = ax.vlines(-np.median(sens_dataIO), 0, 3500, linestyle=":", lw=2.5, color="red")
    l2 = ax.vlines(-np.median(sens_dataIO_low), 0, 3500, linestyle="-.", lw=2.5, color="orange")
    l3 = ax.vlines(-np.median(sens_dataIO_high), 0, 3500, linestyle="--", lw=2.5, color="brown")
    #l3 = ax.vlines(-asimovIO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(-np.median(sens_dataIO)-4, 1000, f"{-np.median(sens_dataIO):.1f}", fontsize=15, color="red")
    ax.text(-np.median(sens_dataIO_low)-4, 2000, f"{-np.median(sens_dataIO_low):.1f}", fontsize=15, color="orange")
    ax.text(-np.median(sens_dataIO_high)-4, 2800, f"{-np.median(sens_dataIO_high):.1f}", fontsize=15, color="brown")
    #ax.text(-asimovIO+2, 1500, f"{-asimovIO:.1f}", fontsize=15, color="black")

    ax.vlines(np.median(sens_dataNO),  0, 3500, linestyle=":", lw=2.5, color="blue")
    ax.vlines(np.median(sens_dataNO_low),  0, 3500, linestyle="-.", lw=2.5, color="green")
    ax.vlines(np.median(sens_dataNO_high),  0, 3500, linestyle="--", lw=2.5, color="darkviolet")
    #ax.vlines(asimovNO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(np.median(sens_dataNO)+2, 1000, f"{np.median(sens_dataNO):.1f}", fontsize=15, color="blue")
    ax.text(np.median(sens_dataNO_low)+2, 2000, f"{np.median(sens_dataNO_low):.1f}", fontsize=15, color="green")
    ax.text(np.median(sens_dataNO_high)+2, 2800, f"{np.median(sens_dataNO_high):.1f}", fontsize=15, color="darkviolet")
    #ax.text(asimovNO+2, 1500, f"{asimovNO:.1f}", fontsize=15, color="black")
    ax.set_ylim(0, 7000)
    
    #lg1 = ax.legend([l1, l3], ["median value of Monte Carlo", "Asimov dataset"], loc="upper right", frameon=True)
    #lg1 = ax.legend([l1, l2], ["high C14", "low C14"], loc="upper right", frameon=True)
    #lg1 = ax.legend([l1, l2, l3], ["unbinned", "0.01 ms binned", "Asimov"], loc="upper right", frameon=True)
    plt.legend(prop={"size":12}, loc="upper left", frameon=True, ncol=3)
    #plt.gca().add_artist(lg1)

    plt.tight_layout()
    plt.savefig(f"./plots/Garching82703_pESeESIBD_10kpc_{Et:.2f}MeV_withC14bkg.png")
    plt.show()



def draw_Tmax():
    fig, ax = plt.subplots(figsize=(9, 4))
    sens_NO_arr, sens_IO_arr = [], []
    Tmax_arr = [10, 15, 20, 25, 30]
    for tmax in Tmax_arr:
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", 0.1, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
        sens_NO_arr.append(np.median(sens_dataNO))
        sens_IO_arr.append(np.median(sens_dataIO))

    ax.plot(Tmax_arr, sens_NO_arr, "v-", ms=8, lw=2, label="Normal ordering")
    ax.plot(Tmax_arr, sens_IO_arr, "^-", ms=8, lw=2, label="Inverted ordering")


    ax.set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=20)
    ax.set_ylabel(r"$\Delta\chi^2$", fontsize=20)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    ax.grid(True, linestyle=":")
    ax.legend(prop={"size":16}, loc="upper right", frameon=True)

    plt.tight_layout()
    #plt.savefig(f"./plots/Garching82703_pESeESIBD+CEvNS_10kpc_noC14bkg.png")
    plt.show()


def draw_Ethr():
    currentDateAndTime = datetime.now()

    fig, ax = plt.subplots(figsize=(9, 5))
    #fig, (ax1, ax) = plt.subplots(2, 1, figsize=(9, 5))

    Et_arr = [0.00, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, ]
    sens_NO_arr, sens_IO_arr = [], []
    for Et in Et_arr:
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
        sens_NO_arr.append(np.median(sens_dataNO))
        sens_IO_arr.append(np.median(sens_dataIO))

    ax.plot(Et_arr, sens_NO_arr, "^-", ms=8, color="blue", lw=2, label="NO: w/o C14")
    ax.plot(Et_arr, sens_IO_arr, "v-", ms=8, color="crimson", lw=2, label="IO: w/o C14")

    Et_arr = np.arange(0, 0.16, 0.01)
    sens_NO_arrl, sens_IO_arrl = [], []
    for Et in Et_arr:
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_C14bkglow_PosiToyData")
        sens_NO_arrl.append(np.median(sens_dataNO))
        sens_IO_arrl.append(np.median(sens_dataIO))

    sens_NO_arrh, sens_IO_arrh = [], []
    for Et in Et_arr:
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_C14bkghigh_PosiToyData")
        sens_NO_arrh.append(np.median(sens_dataNO))
        sens_IO_arrh.append(np.median(sens_dataIO))

    ax.fill_between(Et_arr, sens_NO_arrl, sens_NO_arrh, color='blue', alpha=0.5, label="NO: w/ C14")
    ax.fill_between(Et_arr, sens_IO_arrl, sens_IO_arrh, color="crimson", alpha=0.5, label="IO: w/ C14")


    #_, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBDCEvNS", 0.00, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #ax.plot(0, np.median(sens_dataNO), "v", ms=8, color="darkviolet")
    #ax.plot(0, np.median(sens_dataIO), "^", ms=8, color="darkviolet")
    #_, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", 0.10, 20, 10, end="tot_nodoFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #ax.hlines(np.median(sens_dataNO), 0, 0.2, linestyle="--", lw=2, alpha=0.6, color="green")
    #ax.hlines(np.median(sens_dataIO), 0, 0.2, linestyle="-.", lw=2, alpha=0.6, color="green")


    #ax.semilogx()
    ax.set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=18)
    ax.set_ylabel(r"$\Delta\chi^2_\mathrm{m}$", fontsize=18)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.grid(True, linestyle=":")
    ax.legend(prop={"size":12}, ncol=2, loc="upper right", frameon=True)
    ax.set_ylim(2, 11)

    """
    MO = "NO"
    Nevt = [[] for  i in range(3)]
    Nevt_data = [[] for  i in range(3)]
    NevtErr_data = [[] for i in range(3)]
    #Ethr = np.arange(0.05, 0.16, 0.01)
    Ethr = Et_arr
    for E in Ethr:
        for i, cha in enumerate(["eES", "pES", "IBD"]):
            filename = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{MO}_10kpc_{cha}_{E:.2f}MeV_newshortPDF.root"
            print(filename)
            f = ROOT.TFile(filename, "read")
            h = f.Get("h1")
            bin1, bin2 = h.FindBin(-20), h.FindBin(20)
            N = h.Integral(bin1, bin2, "width")
            Nevt[i].append(N)

            #filename = f"./scale1_poisson/Garching82703_{cha}_binneddata_{MO}_10kpc_thr{E:.2f}MeV_Tmin-20msTmax20ms_binning_new.root"
            #f = up.open(filename)
            #arr = f["binned"]["TbinConts"].array()
            #NperEvt = []
            #for subarr in arr:
            #    NperEvt.append(len(subarr))
            #NperEvt = np.array(NperEvt)
            #Nevt_data[i].append(np.mean(NperEvt)) 
            #NevtErr_data[i].append(np.std(NperEvt))


    ax1.plot(Ethr, Nevt[0], "^-", ms=5, lw=2, fillstyle="none", color="blue", label="eES (NO)")
    ax1.plot(Ethr, Nevt[2], "^-", ms=5, lw=2, fillstyle="none", color="darkviolet", label="IBD (NO)")
    ax1.plot(Ethr, Nevt[1], "s-", ms=5, lw=2, color="red", label="pES")
    MO = "IO"
    Nevt = [[] for  i in range(3)]
    Nevt_data = [[] for  i in range(3)]
    NevtErr_data = [[] for i in range(3)]
    #Ethr = np.arange(0.05, 0.16, 0.01)
    Ethr = Et_arr
    for E in Ethr:
        for i, cha in enumerate(["eES", "pES", "IBD"]):
            filename = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{MO}_10kpc_{cha}_{E:.2f}MeV_newshortPDF.root"
            print(filename)
            f = ROOT.TFile(filename, "read")
            h = f.Get("h1")
            bin1, bin2 = h.FindBin(-20), h.FindBin(20)
            N = h.Integral(bin1, bin2, "width")
            Nevt[i].append(N)
    #ax1.errorbar(Ethr, Nevt_data[0], yerr=NevtErr_data[0], fmt="s-", ms=10, lw=2, fillstyle="none", color="blue", label="eES-data")
    #ax1.errorbar(Ethr, Nevt_data[1], yerr=NevtErr_data[1], fmt="s-", ms=10, lw=2, fillstyle="none", color="red", label="pES-data")
    #ax1.errorbar(Ethr, Nevt_data[2], yerr=NevtErr_data[2], fmt="s-", ms=10, lw=2, fillstyle="none", color="darkviolet", label="IBD-data")

    ax1.plot(Ethr, Nevt[0], "s:", ms=5, lw=2, fillstyle="none", color="blue", label="eES (IO)")
    ax1.plot(Ethr, Nevt[2], "s:", ms=5, lw=2, fillstyle="none", color="darkviolet", label="IBD (IO)")
    ax1.set_ylabel(r"$N_\mathrm{evt}$", fontsize=16)
    ax1.tick_params(axis="y", labelsize=15)
    ax1.grid(True, linestyle=":")
    ax1.legend(prop={"size":12}, loc="lower left", frameon=True, ncol=2)
    ax1.semilogx()
    ax1.semilogy()
    """

    plt.subplots_adjust(hspace=0.0, bottom=0.2)
    plt.savefig(f"../plots/Garching82703_pESeESIBD+CEvNS_10kpc_C14bkg_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}.pdf")
    plt.show()

def draw_models():

    #fig = plt.figure(figsize=(12, 6))
    #spec = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[1, 2])
    #ax0 = fig.add_subplot(spec[0])
    #ax = fig.add_subplot(spec[1])
    
    fig, ax = plt.subplots(figsize=(10, 5))
    colors0 = mpl.cm.gnuplot(np.linspace(0.1,1, 8))
    colors1 = mpl.cm.viridis(np.linspace(0.1,1, 8))
    sens_NO_arr, sens_IO_arr = [], []
    mods = ["81123", "82503", "82703", "84003", "91123", "92503", "92703", "94003"]
    labels = [r"Shen EoS, $11.2M_\odot$", r"Shen EoS, $25M_\odot$", r"Shen EoS, $27M_\odot$", r"Shen EoS, $40M_\odot$", r"LS200 EoS, $11.2M_\odot$", r"LS200 EoS, $25M_\odot$", r"LS200 EoS, $27M_\odot$", r"LS200 EoS, $40M_\odot$", ]
    label1 = [r"$12M_\odot$", r"$14M_\odot$", r"$16M_\odot$", r"$18M_\odot$", r"$20M_\odot$", r"$22M_\odot$", r"$25M_\odot$", r"$26M_\odot$", ]
    tmax_arr = [15, 20, 25, 30]
    for i, mod in enumerate(mods):
        m_sens_dataNO, m_sens_dataIO = [], []
        for tmax in tmax_arr:
            _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching"+mod, "pESeESIBD", 0.10, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
            m_sens_dataNO.append(np.median(sens_dataNO))
            m_sens_dataIO.append(np.median(sens_dataIO))

        ax.plot(tmax_arr, m_sens_dataNO, "^", color=colors0[i], ms=8, lw=2, label=labels[i])
        ax.plot(tmax_arr, m_sens_dataIO, "^", color=colors0[i], ms=8, lw=2, fillstyle="none")


    #mods = ["82503", "82703", "84003", "91123", "92503", "92703", "94003"]
    #stdNO_Garching, stdIO_Garching = [], []
    #onetNO, onetIO = [], []
    #for mod in mods:
    #    _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching"+mod, "pESeES", 0.10, 10, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #    onetNO.append(np.median(sens_dataNO))
    #    onetIO.append(np.median(sens_dataIO))
    #onetNO = np.array(onetNO)
    #onetIO = np.array(onetIO)
    #stdNO_Garching.append(np.std(onetNO))
    #stdIO_Garching.append(np.std(onetIO))
    #mods = ["81123", "82503", "82703", "84003", "91123", "92503", "92703", "94003"]
    #for tmax in tmax_arr:
    #    onetNO, onetIO = [], []
    #    for mod in mods:
    #        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching"+mod, "pESeESIBD", 0.10, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #        onetNO.append(np.median(sens_dataNO))
    #        onetIO.append(np.median(sens_dataIO))
    #    onetNO = np.array(onetNO)
    #    onetIO = np.array(onetIO)
    #    stdNO_Garching.append(np.std(onetNO))
    #    stdIO_Garching.append(np.std(onetIO))
    #ax0.plot(np.arange(10, 31, 5), stdNO_Garching, "^-", ms=8, lw=2, label="Garching, NO")
    #ax0.plot(np.arange(10, 31, 5), stdIO_Garching, "^-", ms=8, lw=2, label="Garching, IO", fillstyle="none")

    t_bur = np.array([45, 50, 55, 60])
    for i , mod in enumerate([12, 14, 16, 18, 20, 22, 25, 26]):
        m_sens_dataNO, m_sens_dataIO = [], []
        for tmax in t_bur:
            _, _, sens_dataNO, _, _, sens_dataIO            = read("Burrows"+str(mod), "pESeESIBD", 0.10, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
            m_sens_dataNO.append(np.median(sens_dataNO))
            m_sens_dataIO.append(np.median(sens_dataIO))

        ax.plot(t_bur-30, m_sens_dataNO, "s", color=colors1[i], ms=8, lw=2, label=label1[i])
        ax.plot(t_bur-30, m_sens_dataIO, "s", color=colors1[i], ms=8, lw=2, fillstyle="none")

    #_, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBDCEvNS", 0.00, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #ax.plot(0, np.median(sens_dataNO), "v", ms=8, color="red")
    #ax.plot(0, np.median(sens_dataIO), "^", ms=8, color="red")

    #stdNO_Bur, stdIO_Bur = [], []
    #for tmax in [40]:
    #    onetNO, onetIO = [], []
    #    for mod in [12, 14, 16, 18, 20, 22, 25, 26]:
    #        _, _, sens_dataNO, _, _, sens_dataIO            = read("Burrows"+str(mod), "pESeES", 0.10, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #        onetNO.append(np.median(sens_dataNO))
    #        onetIO.append(np.median(sens_dataIO))
    #    onetNO = np.array(onetNO)
    #    onetIO = np.array(onetIO)
    #    stdNO_Bur.append(np.std(onetNO))
    #    stdIO_Bur.append(np.std(onetIO))
    #for tmax in t_bur:
    #    onetNO, onetIO = [], []
    #    for mod in [12, 14, 16, 18, 20, 22, 25, 26]:
    #        _, _, sens_dataNO, _, _, sens_dataIO            = read("Burrows"+str(mod), "pESeESIBD", 0.10, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #        onetNO.append(np.median(sens_dataNO))
    #        onetIO.append(np.median(sens_dataIO))
    #    onetNO = np.array(onetNO)
    #    onetIO = np.array(onetIO)
    #    stdNO_Bur.append(np.std(onetNO))
    #    stdIO_Bur.append(np.std(onetIO))
    #ax0.plot(np.arange(10, 31, 5), stdNO_Bur, "s-", ms=8, lw=2, label="Burrows, NO")
    #ax0.plot(np.arange(10, 31, 5), stdIO_Bur, "s-", ms=8, lw=2, label="Burrows, IO", fillstyle="none")


    sens_NO_arr, sens_IO_arr = [], []
    sens_NO_arr.append(0)
    sens_IO_arr.append(0)
    mods = ["82503", "82703", "84003", "91123", "92503", "92703", "94003"]
    for i, mod in enumerate(mods):
        m_sens_dataNO, m_sens_dataIO = [], []
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching"+mod, "pESeES", 0.10, 10, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
        m_sens_dataNO.append(np.median(sens_dataNO))
        m_sens_dataIO.append(np.median(sens_dataIO))
            
        ax.plot([10], m_sens_dataNO, "^", color=colors0[i], ms=8, lw=2)
        ax.plot([10], m_sens_dataIO, "^", color=colors0[i], ms=8, lw=2, fillstyle="none")

    for i , mod in enumerate([12, 14, 16, 18, 20, 22, 25, 26]):
        m_sens_dataNO, m_sens_dataIO = [], []
        for tmax in [40]:
            _, _, sens_dataNO, _, _, sens_dataIO            = read("Burrows"+str(mod), "pESeES", 0.10, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
            m_sens_dataNO.append(np.median(sens_dataNO))
            m_sens_dataIO.append(np.median(sens_dataIO))
        ax.plot([10], m_sens_dataNO, "s", color=colors1[i], ms=8, lw=2)
        ax.plot([10], m_sens_dataIO, "s", color=colors1[i], ms=8, lw=2, fillstyle="none")



    currentDateAndTime = datetime.now()


    #ax.semilogx()
    #ax0.grid(True)
    #ax0.set_ylabel(r"$\sigma(\Delta\chi^2_\mathrm{m})$", fontsize=18)
    #ax0.tick_params(axis="x", labelsize=0)
    #ax0.tick_params(axis="y", labelsize=16)
    #ax0.set_xlim(9, 31)
    ax.set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=18)
    ax.set_ylabel(r"$\Delta\chi^2_\mathrm{m}$", fontsize=18)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.grid(True, linestyle=":")
    ax.legend(prop={"size":11}, loc="lower right", frameon=True, ncol=4)
    ax.set_xlim(9, 31)
    #ax0.legend(prop={"size":13}, frameon=True, ncol=2)

    #ax.text(9.5, 4.5, "Normal ordering", fontsize=15)
    #ax.text(9.5, 6.0, "Inveted ordering", fontsize=15)

    plt.tight_layout()
    #plt.subplots_adjust(hspace=0)
    plt.savefig(f"../plots/Garchings+Burrows_pESeESIBD_202530ms_10kpc_noC14bkg_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}.pdf")
    plt.show()


def drawStat():
    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO, NeES_NO, NpES_NO, NIBD_NO, NeES_IO, NpES_IO, NIBD_IO = read2("../results/Garching82703_10kpc_NO_pESeESIBD_useMassFalse_nuMass0.0eV_0.60MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_THEIA100.csv", "../results/Garching82703_10kpc_IO_pESeESIBD_useMassFalse_nuMass0.0eV_0.60MeV_fitTmin-0.020sfitTmax0.020s_data1D_tot_unbinnedData_unbinnedNLL_noC14bkg_PosiToyDataTobs1D_2023313_THEIA100.csv")

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(NeES_NO, bins=20, alpha=0.5, range=(0, 20), label=f"eES: {np.mean(NeES_NO):.1f}")
    ax[0].hist(NIBD_NO, bins=50, alpha=0.5, range=(0, 150), label=f"IBD: {np.mean(NIBD_NO):.1f}")
    ax[0].hist(NpES_NO, bins=50, alpha=0.5, range=(0, 50), label=f"pES: {np.mean(NpES_NO):.1f}")

    ax[1].hist(NeES_IO, bins=20, alpha=0.5, range=(0, 20), label=f"eES: {np.mean(NeES_IO):.1f}")
    ax[1].hist(NIBD_IO, bins=50, alpha=0.5, range=(0, 150), label=f"IBD: {np.mean(NIBD_IO):.1f}")
    ax[1].hist(NpES_IO, bins=50, alpha=0.5, range=(0, 50), label=f"pES: {np.mean(NpES_IO):.1f}")

    ax[0].set_xlabel("event number", fontsize=17)
    ax[0].tick_params(axis="x", labelsize=15)
    ax[0].legend(prop={"size":14})
    ax[1].set_xlabel("event number", fontsize=17)
    ax[1].tick_params(axis="x", labelsize=15)
    ax[1].legend(prop={"size":14})

    plt.tight_layout()
    plt.show()
    


def drawNevent():
    MO = "IO"
    Nevt = [[] for  i in range(3)]
    Nevt_data = [[] for  i in range(3)]
    NevtErr_data = [[] for i in range(3)]
    #Ethr = np.arange(0.05, 0.16, 0.01)
    Ethr = [0.10]
    for E in Ethr:
        for i, cha in enumerate(["eES", "pES", "IBD"]):
            filename = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{MO}_10kpc_{cha}_{E:.2f}MeV_newshortPDF.root"
            print(filename)
            f = ROOT.TFile(filename, "read")
            h = f.Get("h1")
            bin1, bin2 = h.FindBin(-20), h.FindBin(20)
            N = h.Integral(bin1, bin2, "width")
            Nevt[i].append(N)

            filename = f"./scale1_poisson/Garching82703_{cha}_binneddata_{MO}_10kpc_thr{E:.2f}MeV_Tmin-20msTmax20ms_binning_new.root"
            f = up.open(filename)
            arr = f["binned"]["TbinConts"].array()
            NperEvt = []
            for subarr in arr:
                NperEvt.append(len(subarr))
            NperEvt = np.array(NperEvt)
            Nevt_data[i].append(np.mean(NperEvt)) 
            NevtErr_data[i].append(np.std(NperEvt))


    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(Ethr, Nevt[0], "o-", ms=8, lw=2, color="blue", label="eES-pdf")
    ax.plot(Ethr, Nevt[1], "o-", ms=8, lw=2, color="red", label="pES-pdf")
    ax.plot(Ethr, Nevt[2], "o-", ms=8, lw=2, color="darkviolet", label="IBD-pdf")
    ax.errorbar(Ethr, Nevt_data[0], yerr=NevtErr_data[0], fmt="s-", ms=10, lw=2, fillstyle="none", color="blue", label="eES-data")
    ax.errorbar(Ethr, Nevt_data[1], yerr=NevtErr_data[1], fmt="s-", ms=10, lw=2, fillstyle="none", color="red", label="pES-data")
    ax.errorbar(Ethr, Nevt_data[2], yerr=NevtErr_data[2], fmt="s-", ms=10, lw=2, fillstyle="none", color="darkviolet", label="IBD-data")
    ax.set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax.set_ylabel(r"expected event number", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.grid(True, linestyle=":")
    ax.legend(prop={"size":12}, loc="upper right", frameon=True, ncol=3)

    plt.tight_layout()
    #plt.savefig("./plots/Garching82703_10kpc_diffEthr_IO_Nevt.pdf")
    plt.show()


def draw_distance():

    fig, (ax, ax1) = plt.subplots(1, 2, figsize=(12, 5))
    #colors0 = mpl.cm.gnuplot(np.linspace(0.1,1, 8))
    sens_NO_arr, sens_IO_arr = [], []
    dist = [5, 6, 7, 8, 9, 10]
    ndist = len(dist)
    Ethrs = [0.10, 0.11, 0.12, 0.13, 0.14, 0.15]
    nEthrs = len(Ethrs)
    m_sens_dataNO = np.zeros((nEthrs, ndist))
    m_sens_dataIO = np.zeros((nEthrs, ndist))
    for i, d in enumerate(dist):
        for j, E in enumerate(Ethrs):
            _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", E, 20, d, end="tot_doFit_binnedData_unbinnedNLL_C14bkglow_PosiToyData")
            m_sens_dataNO[nEthrs-1-j, i] = (np.median(sens_dataNO))
            m_sens_dataIO[nEthrs-1-j, i] = (np.median(sens_dataIO))
        
    #ax.plot(dist, m_sens_dataNO, "^", ms=8, lw=2)
    #ax.plot(dist, m_sens_dataIO, "v", ms=8, lw=2, fillstyle="none")

    X = np.arange(5, 11, 1)
    Y = np.arange(0.10, 0.16, 0.01)

    im = ax.imshow(m_sens_dataNO,   extent=[4.5, 10.5, 0.095, 0.155], aspect="auto")
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(r"$\Delta \chi^2_\mathrm{m}$", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    im1 = ax1.imshow(m_sens_dataIO, extent=[4.5, 10.5, 0.095, 0.155], aspect="auto")
    cb1 = plt.colorbar(im1, ax=ax1)
    cb1.ax.tick_params(labelsize=14)
    cs = ax.contour(X, Y, m_sens_dataNO, levels=[9, 16,  25], colors=[ "red", "darkorange", "darkviolet"], linewidth=2)
    cs1 = ax1.contour(X, Y, m_sens_dataIO ,levels=[9, 16, 25],colors=[ "red", "darkorange", "darkviolet"], linewidth=2)

    ax.text(5.0, 0.15, r"$\Delta\chi^2_\mathrm{m}$ = 25",  color="darkviolet", fontsize=14)
    ax.text(6.5, 0.13, r"$\Delta\chi^2_\mathrm{m}$ = 16",  color="darkorange", fontsize=14)
    ax.text(8.5, 0.11, r"$\Delta\chi^2_\mathrm{m}$ = 9",   color="red", fontsize=14)
    ax1.text(5.0, 0.15, r"$\Delta\chi^2_\mathrm{m}$ = 25",  color="darkviolet", fontsize=14)
    ax1.text(7.0, 0.13, r"$\Delta\chi^2_\mathrm{m}$ = 16",  color="darkorange", fontsize=14)
    ax1.text(9.2, 0.11, r"$\Delta\chi^2_\mathrm{m}$ = 9",   color="red", fontsize=14)

    cb1.set_label(r"$\Delta \chi^2_\mathrm{m}$", fontsize=14)

    #ax.semilogx()
    ax.set_xlabel("Distance [kpc]", fontsize=16)
    ax.set_ylabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax.set_xticks(np.arange(5, 11 ,1), ["5", "6", "7", "8", "9", "10"])
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax1.set_xticks(np.arange(5, 11 ,1), ["5", "6", "7", "8", "9", "10"])
    ax1.set_xlabel("Distance [kpc]", fontsize=16)
    ax1.set_ylabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax1.tick_params(axis="x", labelsize=15)
    ax1.tick_params(axis="y", labelsize=15)


    currentDateAndTime = datetime.now()
    plt.tight_layout()
    plt.savefig(f"../plots/Garching82703_pESeESIBD_20ms_distances_C14bkglow_{currentDateAndTime.year}{currentDateAndTime.month}{currentDateAndTime.day}.pdf")
    plt.show()


if __name__ == "__main__":
    #draw()
    #drawHist2d()
    #drawNLL()
    #scatter()
    #drawFitT()
    #drawStat()
    #drawC14()
    #draw_Ethr()
    #drawNevent()
    #draw_Tmax()
    draw_models()
    #draw_distance()




