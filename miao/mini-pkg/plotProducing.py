import numpy as np
import matplotlib.pyplot as plt
import SNGarchingIntegFcn as gar
from SNnumGarchingSrc import SNnumGarchingSrc

imode = 82500
gar.readFluxGraph(imode)
nuTypeName = [r"$\nu_e$", r"$\nu_{\bar{e}}$", r"$\nu_x$", r"$\nu_{\bar{x}}$"]

def Lum_vs_time_allType():
    time_nue_burst, lum_nue_burst = [], []
    time_nue_acc, lum_nue_acc = [], []
    time_nue_cool, lum_nue_cool = [], []
    time_nuebar_burst, lum_nuebar_burst = [], []
    time_nuebar_acc, lum_nuebar_acc = [], []
    time_nuebar_cool, lum_nuebar_cool = [], []
    time_nux_burst, lum_nux_burst = [], []
    time_nux_acc, lum_nux_acc = [], []
    time_nux_cool, lum_nux_cool = [], []

    for i in range(gar.grLuminosity[0].GetN()):
        t = gar.grLuminosity[0].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nue_burst.append(gar.grLuminosity[0].GetPointX(i))
            lum_nue_burst.append(gar.grLuminosity[0].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nue_acc.append(gar.grLuminosity[0].GetPointX(i))
            lum_nue_acc.append(gar.grLuminosity[0].GetPointY(i))
        if 0.8 < t < 8:
            time_nue_cool.append(gar.grLuminosity[0].GetPointX(i))
            lum_nue_cool.append(gar.grLuminosity[0].GetPointY(i))

    for i in range(gar.grLuminosity[1].GetN()):
        t = gar.grLuminosity[1].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nuebar_burst.append(gar.grLuminosity[1].GetPointX(i))
            lum_nuebar_burst.append(gar.grLuminosity[1].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nuebar_acc.append(gar.grLuminosity[1].GetPointX(i))
            lum_nuebar_acc.append(gar.grLuminosity[1].GetPointY(i))
        if 0.8 < t < 8:
            time_nuebar_cool.append(gar.grLuminosity[1].GetPointX(i))
            lum_nuebar_cool.append(gar.grLuminosity[1].GetPointY(i))


    for i in range(gar.grLuminosity[2].GetN()):
        t = gar.grLuminosity[2].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nux_burst.append(gar.grLuminosity[2].GetPointX(i))
            lum_nux_burst.append(gar.grLuminosity[2].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nux_acc.append(gar.grLuminosity[2].GetPointX(i))
            lum_nux_acc.append(gar.grLuminosity[2].GetPointY(i))
        if 0.8 < t < 8:
            time_nux_cool.append(gar.grLuminosity[2].GetPointX(i))
            lum_nux_cool.append(gar.grLuminosity[2].GetPointY(i))

    plt.figure(0)
    plt.plot(time_nue_burst, lum_nue_burst, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_burst, lum_nuebar_burst, "-", label=nuTypeName[1])
    plt.plot(time_nux_burst, lum_nux_burst, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"Luminosity $10^{51}$ erg/s")
    plt.savefig("./outputs/Lum_vs_time_allType_burst_Garching%d.pdf"%imode)

    plt.figure(1)
    plt.plot(time_nue_acc, lum_nue_acc, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_acc, lum_nuebar_acc, "-", label=nuTypeName[1])
    plt.plot(time_nux_acc, lum_nux_acc, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"Luminosity $10^{51}$ erg/s")
    plt.savefig("./outputs/Lum_vs_time_allType_accretion_Garching%d.pdf"%imode)

    plt.figure(2)
    plt.plot(time_nue_cool, lum_nue_cool, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_cool, lum_nuebar_cool, "-", label=nuTypeName[1])
    plt.plot(time_nux_cool, lum_nux_cool, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"Luminosity $10^{51}$ erg/s")
    plt.savefig("./outputs/Lum_vs_time_allType_cooling_Garching%d.pdf" %imode)

    print("Lum_vs_time_allType.pdf has been created !!!")


def AveE_vs_time_allType():

    time_nue_burst, aveE_nue_burst = [], []
    time_nue_acc, aveE_nue_acc = [], []
    time_nue_cool, aveE_nue_cool = [], []
    time_nuebar_burst, aveE_nuebar_burst = [], []
    time_nuebar_acc, aveE_nuebar_acc = [], []
    time_nuebar_cool, aveE_nuebar_cool = [], []
    time_nux_burst, aveE_nux_burst = [], []
    time_nux_acc, aveE_nux_acc = [], []
    time_nux_cool, aveE_nux_cool = [], []

    for i in range(gar.grAverageE[0].GetN()):
        t = gar.grAverageE[0].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nue_burst.append(gar.grAverageE[0].GetPointX(i))
            aveE_nue_burst.append(gar.grAverageE[0].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nue_acc.append(gar.grAverageE[0].GetPointX(i))
            aveE_nue_acc.append(gar.grAverageE[0].GetPointY(i))
        if 0.8 < t < 8:
            time_nue_cool.append(gar.grAverageE[0].GetPointX(i))
            aveE_nue_cool.append(gar.grAverageE[0].GetPointY(i))

    for i in range(gar.grAverageE[1].GetN()):
        t = gar.grAverageE[1].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nuebar_burst.append(gar.grAverageE[1].GetPointX(i))
            aveE_nuebar_burst.append(gar.grAverageE[1].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nuebar_acc.append(gar.grAverageE[1].GetPointX(i))
            aveE_nuebar_acc.append(gar.grAverageE[1].GetPointY(i))
        if 0.8 < t < 8:
            time_nuebar_cool.append(gar.grAverageE[1].GetPointX(i))
            aveE_nuebar_cool.append(gar.grAverageE[1].GetPointY(i))


    for i in range(gar.grAverageE[2].GetN()):
        t = gar.grAverageE[2].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nux_burst.append(gar.grAverageE[2].GetPointX(i))
            aveE_nux_burst.append(gar.grAverageE[2].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nux_acc.append(gar.grAverageE[2].GetPointX(i))
            aveE_nux_acc.append(gar.grAverageE[2].GetPointY(i))
        if 0.8 < t < 8:
            time_nux_cool.append(gar.grAverageE[2].GetPointX(i))
            aveE_nux_cool.append(gar.grAverageE[2].GetPointY(i))

    plt.figure(0)
    plt.plot(time_nue_burst, aveE_nue_burst, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_burst, aveE_nuebar_burst, "-", label=nuTypeName[1])
    plt.plot(time_nux_burst, aveE_nux_burst, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("AverageE/MeV")
    plt.savefig("./outputs/AveE_vs_time_allType_burst_Garching%d.pdf" %imode)

    plt.figure(1)
    plt.plot(time_nue_acc, aveE_nue_acc, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_acc, aveE_nuebar_acc, "-", label=nuTypeName[1])
    plt.plot(time_nux_acc, aveE_nux_acc, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("AverageE/MeV")
    plt.savefig("./outputs/AveE_vs_time_allType_accretion_Garching%d.pdf"%imode)

    plt.figure(2)
    plt.plot(time_nue_cool, aveE_nue_cool, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_cool, aveE_nuebar_cool, "-", label=nuTypeName[1])
    plt.plot(time_nux_cool, aveE_nux_cool, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("AverageE/MeV")
    plt.savefig("./outputs/AveE_vs_time_allType_cooling_Garching%d.pdf"%imode)

    print("AveE_vs_time_allType.pdf has been created !!!")





