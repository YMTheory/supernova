import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import units as u

from CEvNS_XS import CEvNS_XS
from IBD_XS import IBD_XS
from NuE_XS import NuE_XS
from NuP_XS import NuP_XS
from SNNuloader import SNNuloader
from toyDetector import toyDetector

from visible_spectrum import visible_spectrum

from tqdm import  tqdm

import sys

def drawMO_ufunc(cha, SNmass=25, EoS="LS220"):
    spec = visible_spectrum(cha, SNmass=SNmass, EoS=EoS)
    tarr = np.arange(-0.03, 0.03, 0.001)

    fig, ax = plt.subplots(figsize=(6, 4))

    Evismin, Evismax = 0.10, 80
    no = np.zeros(len(tarr))
    io = np.zeros(len(tarr))
    for i in tqdm(range(len(tarr))):
        t = tarr[i]
        nevt_no = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "no")
        no[i] = nevt_no
        nevt_io = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "io")
        io[i] = nevt_io

    ax.plot(tarr*1000, no/1000., "-",  lw=2, label="NO")
    ax.plot(tarr*1000, io/1000., "-",  lw=2, label="IO")

    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("counts per ms", fontsize=14)
    ax.legend(prop={"size":12})
    ax.tick_params(axis="both", labelsize=13)
    ax.grid(True, linestyle=":")
    plt.tight_layout()
    plt.savefig(f"./plots/timespec_{cha}_{EoS}{SNmass}.pdf")
    plt.show()

    return tarr, no, io



def drawEthr_ufunc(cha, SNmass=25, EoS="LS220"):
    spec = visible_spectrum(cha, SNmass=SNmass, EoS=EoS)
    tarr = np.arange(-0.03, 0.03, 0.001)

    fig, ax = plt.subplots(figsize=(6, 4))

    Evismin, Evismax = 0.0, 10
    no = np.zeros(len(tarr))
    io = np.zeros(len(tarr))
    for i in tqdm(range(len(tarr))):
        t = tarr[i]
        nevt_no = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "no")
        no[i] = nevt_no
        nevt_io = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "io")
        io[i] = nevt_io

    ax.plot(tarr*1000, no/1000., "-",  lw=2, label=r"E$_\mathrm{thr}$=0.00 MeV")

    Evismin, Evismax = 0.10, 10
    no = np.zeros(len(tarr))
    io = np.zeros(len(tarr))
    for i in tqdm(range(len(tarr))):
        t = tarr[i]
        nevt_no = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "no")
        no[i] = nevt_no
        nevt_io = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "io")
        io[i] = nevt_io

    ax.plot(tarr*1000, no/1000., "-",  lw=2, label=r"E$_\mathrm{thr}$=0.10 MeV")


    Evismin, Evismax = 0.15, 10
    no = np.zeros(len(tarr))
    io = np.zeros(len(tarr))
    for i in tqdm(range(len(tarr))):
        t = tarr[i]
        nevt_no = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "no")
        no[i] = nevt_no
        nevt_io = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "io")
        io[i] = nevt_io

    ax.plot(tarr*1000, no/1000., "-",  lw=2, label=r"E$_\mathrm{thr}$=0.15 MeV")


    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("counts per ms", fontsize=14)
    if cha == "CEvNS":
        ax.ticklabel_format(style='sci', axis='y')
        ax.semilogy()
    ax.legend(prop={"size":12})
    ax.grid(True, linestyle=":")
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()
    plt.savefig(f"./plots/timespec_{cha}_{EoS}{SNmass}.pdf")
    plt.show()

    return tarr, no, io



def test(Emin, Emax, cha, SNmass, EoS):

    spec = visible_spectrum(cha, SNmass=SNmass, EoS=EoS)
    tarr = np.arange(-0.03, 0.03, 0.001)

    fig, ax = plt.subplots(figsize=(6, 4))

    no = np.zeros((len(Emax), len(tarr)))
    for j in range(len(Emax)):
        E0 = Emin[j]
        E = Emax[j]
        for i in tqdm(range(len(tarr))):
            t = tarr[i]
            no[j, i] = spec.getVisibleEventAtTEvisIntegral(t, E0, E, "no")


        ax.plot(tarr*1000, no[j]/1000., "-",  lw=2, label=r"E$_\mathrm{max}$"+f"={E} MeV")

    np.savetxt("./data/NuP_check.txt", no)
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("counts per ms", fontsize=14)
    ax.legend(prop={"size":12})
    ax.grid(True, linestyle=":")
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()
    plt.savefig(f"./plots/timespec_{cha}_testInt.pdf")
    plt.show()


def output(cha, SNmass, EoS, Evismin, Evismax, tmin, tmax, tstep):
    spec = visible_spectrum(cha, SNmass=SNmass, EoS=EoS)
    tarr = np.arange(tmin, tmax, tstep)
    arr = np.zeros((len(tarr), 2))
    for j in tqdm(range(len(tarr))):
        t = tarr[j]
        arr[j, 0] = t
        arr[j, 1] = spec.getVisibleEventAtTEvisIntegral(t, Evismin, Evismax, "no") 
        print(j, arr[j, 0], arr[j, 1])

    np.savetxt(f"/junofs/users/miaoyu/supernova/simulation/data/NuP_{Evismin}to{Evismax}MeV_{tmin}to{tmax}s.txt", arr)




if __name__ == '__main__':

    #drawMO_ufunc("IBD", SNmass=25, EoS="LS220")
    drawEthr_ufunc("CEvNS", SNmass=25, EoS="LS220")
    #test([0.1, 0.1, 0.1], [3, 4, 5], "NuP", 25, "LS220")
    #output( "CEvNS", 25, "LS220")

    #tmin = float(sys.argv[1])
    #tmax = float(sys.argv[2])
    #Emax = float(sys.argv[3])
    #output("NuP", 25, "LS220", 0.1, Emax, tmin, tmax, 0.001)
    







