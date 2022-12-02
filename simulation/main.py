import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import units as u

from CEvNS_XS import CEvNS_XS
from IBD_XS import IBD_XS
from NuE_XS import NuE_XS
from SNNuloader import SNNuloader
from toyDetector import toyDetector

from visible_spectrum import visible_spectrum

from tqdm import  tqdm

def draw_EvT2D_atEarth(model, flavorID):
    model.load()
    fnu_noosc  = np.zeros((60, 100))
    fnu_no     = np.zeros((60, 100))
    fnu_io     = np.zeros((60, 100))
    tmin, tmax, stept = -0.02, 0.08, 0.001
    Numt = int((tmax - tmin) / stept)
    Evmin, Evmax, stepEv = 0, 60, 1
    NumEv = int((Evmax - Evmin)/stepEv)
    for i in range(Numt):
        t = tmin + (i+0.5) * stept
        for j in range(NumEv):
            Ev = Evmin + (j + 0.5) * stepEv
            fnu_noosc[NumEv-1-j, i] = model.getEventAtEarthAtTime(t, Ev, "noosc", flavorID)
            fnu_no[NumEv-1-j, i] = model.getEventAtEarthAtTime(t, Ev, "no", flavorID)
            fnu_io[NumEv-1-j, i] = model.getEventAtEarthAtTime(t, Ev, "io", flavorID)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    im0 = ax[0].imshow(fnu_noosc, extent=[-20, 80, 0, 60], aspect="auto")
    plt.colorbar(im0, ax=ax[0], shrink=0.6)
    im1 = ax[1].imshow(fnu_no, extent=[-20, 80, 0, 60], aspect="auto")
    plt.colorbar(im1, ax=ax[1], shrink=0.6)
    im2 = ax[2].imshow(fnu_io, extent=[-20, 80, 0, 60], aspect="auto")
    plt.colorbar(im2, ax=ax[2], shrink=0.6)
    ax[0].set_xlabel("post-bounce time [ms]", fontsize=14)
    ax[0].set_ylabel(r"$E_\nu$ [MeV]", fontsize=14)
    ax[0].tick_params(axis="both", labelsize=13)
    ax[0].set_title("No osc", fontsize=15)
    ax[1].set_xlabel("post-bounce time [ms]", fontsize=14)
    ax[1].set_ylabel(r"$E_\nu$ [MeV]", fontsize=14)
    ax[1].tick_params(axis="both", labelsize=13)
    ax[1].set_title("NO", fontsize=15)
    ax[2].set_xlabel("post-bounce time [ms]", fontsize=14)
    ax[2].set_ylabel(r"$E_\nu$ [MeV]", fontsize=14)
    ax[2].tick_params(axis="both", labelsize=13)
    ax[2].set_title("IO", fontsize=15)

    plt.tight_layout()
    plt.savefig(f"ET2D_flavor{flavorID}_10kpc.pdf")
    plt.show()


def calc_visibleSpectra(model, xsec, det, Evismin=0, Evismax=0.1, chaname="IBD"):
    spec = visible_spectrum(chaname, model, xsec, det)

    t = np.arange(-0.02, 0.04, 0.001)
    Numt = len(t)
    #Evismin, Evismax, stepEvis = 0, 0.1, 0.005
    stepEvis = 0.005
    NumEvis = int((Evismax - Evismin) / stepEvis)

    no, io = np.zeros(Numt), np.zeros(Numt)
    for i in tqdm(range(Numt)):
        ti = t[i]
        valno = spec.getVisibleEventAtEvisAtT(ti, Evismin, Evismax, "no")
        valio = spec.getVisibleEventAtEvisAtT(ti, Evismin, Evismax, "io")
        no[i] = valno
        io[i] = valio
    return t, no, io



def calc_visibleSpectra2D(model, xsec, det, Evismin=0, Evismax=0.1, chaname="IBD"):
    solM = model.mass
    eos = model.EoS
    spec = visible_spectrum(chaname, model, xsec, det)

    t = np.arange(-0.02, 0.04, 0.001)
    #Evismin, Evismax, stepEvis = 0, 0.2, 0.005
    stepEvis = 0.005
    NumEvis = int((Evismax - Evismin) / stepEvis)
    
    no = np.zeros(len(t))
    io = np.zeros(len(t))
    no2d = np.zeros((NumEvis, len(t)))
    io2d = np.zeros((NumEvis, len(t)))
    #for i, ti in tqdm(enumerate(t)):
    for i in tqdm(range(len(t))):
        ti = t[i]
        valno, valio = 0, 0
        for j in tqdm(range(NumEvis)):
            Emin = Evismin + (j) * stepEvis
            Emax = Evismin + (j + 1) * stepEvis
            tmpno = spec.getVisibleEventAtEvisAtT(ti, Emin, Emax, "no")
            tmpio = spec.getVisibleEventAtEvisAtT(ti, Emin, Emax, "io")
            #print(ti, Emin, Emax, tmpno, tmpio)
            no2d[NumEvis-1-j, i] = tmpno
            io2d[NumEvis-1-j, i] = tmpio
            valno += tmpno
            valio += tmpio
        no[i] = valno
        io[i] = valio
    
    return t, no, io, no2d, io2d

def drawIBD():
    det = toyDetector(20000, 0.12, 0.88)
    xsec = IBD_XS()
    dist = 10
    model0 = SNNuloader(11.2, "Accretion", "Shen", dist)
    model1 = SNNuloader(25, "Accretion", "Shen", dist)
    model2 = SNNuloader(27, "Accretion", "Shen", dist)
    model3 = SNNuloader(40, "Accretion", "Shen", dist)
    model4 = SNNuloader(11.2, "Accretion", "LS220", dist)
    model5 = SNNuloader(25, "Accretion", "LS220", dist)
    model6 = SNNuloader(27, "Accretion", "LS220", dist)
    model7 = SNNuloader(40, "Accretion", "LS220", dist)
    #draw_EvT2D_atEarth(model, 3)
    t0, no0, io0, _, _ = calc_visibleSpectra2D(model0, xsec, det)
    t1, no1, io1, _, _ = calc_visibleSpectra2D(model1, xsec, det)
    t2, no2, io2, _, _ = calc_visibleSpectra2D(model2, xsec, det)
    t3, no3, io3, _, _ = calc_visibleSpectra2D(model3, xsec, det)
    t4, no4, io4, _, _ = calc_visibleSpectra2D(model4, xsec, det)
    t5, no5, io5, _, _ = calc_visibleSpectra2D(model5, xsec, det)
    t6, no6, io6, _, _ = calc_visibleSpectra2D(model6, xsec, det)
    t7, no7, io7, _, _ = calc_visibleSpectra2D(model7, xsec, det)

    colors0 = mpl.cm.gnuplot(np.linspace(0.1,1, 8))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t0*1000, no0/1000, "-", color=colors0[0], lw=2, label="11.2solM, Shen")
    ax.plot(t0*1000, io0/1000, ":", color=colors0[0], lw=2, )
    ax.plot(t1*1000, no1/1000, "-", color=colors0[1], lw=2, label="25solM, Shen")
    ax.plot(t1*1000, io1/1000, ":", color=colors0[1], lw=2, )
    ax.plot(t2*1000, no2/1000, "-", color=colors0[2], lw=2, label="27solM, Shen")
    ax.plot(t2*1000, io2/1000, ":", color=colors0[2], lw=2, )
    ax.plot(t3*1000, no3/1000, "-", color=colors0[3], lw=2, label="40solM, Shen")
    ax.plot(t3*1000, io3/1000, ":", color=colors0[3], lw=2, )
    ax.plot(t4*1000, no4/1000, "-", color=colors0[4], lw=2, label="11.2solM, LS220")
    ax.plot(t4*1000, io4/1000, ":", color=colors0[4], lw=2, )
    ax.plot(t5*1000, no5/1000, "-", color=colors0[5], lw=2, label="25solM, LS220")
    ax.plot(t5*1000, io5/1000, ":", color=colors0[5], lw=2, )
    ax.plot(t6*1000, no6/1000, "-", color=colors0[6], lw=2, label="27solM, LS220")
    ax.plot(t6*1000, io6/1000, ":", color=colors0[6], lw=2, )
    ax.plot(t7*1000, no7/1000, "-", color=colors0[7], lw=2, label="40solM, LS220")
    ax.plot(t7*1000, io7/1000, ":", color=colors0[7], lw=2, )
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("counts per ms", fontsize=14)
    ax.legend(prop={"size":12})
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()
    plt.savefig(f"timespec_ibd_xxj.pdf")
    plt.show()


def drawCEvNS():
    det = toyDetector(20000, 0.12, 0.88)
    Na = 6.022e23
    gtoeV = 5.61e32
    mC12 = 12 / Na * gtoeV
    xsec = CEvNS_XS(mC12, 6, 6)
    dist = 10
    model0 = SNNuloader(25, "Accretion", "LS220", dist)

    #t, no, io = calc_visibleSpectra(model0, xsec, det, Evismin=0.1, Evismax=0.5, chaname="CEvNS")
    t, no, io, _, _ = calc_visibleSpectra2D(model0, xsec, det, Evismin=0.10, Evismax=0.5, chaname="CEvNS")
    np.savetxt("CEvNS1d.txt", no, fmt="%.3e", delimiter=" ")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t*1000, no/1000, "-", lw=2, label="11.2solM, Shen")
    ax.plot(t*1000, io/1000, ":", lw=2, )
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("counts per ms", fontsize=14)
    ax.legend(prop={"size":12})
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()
    plt.savefig(f"timespec_CEvNS.pdf")
    plt.show()


def drawCEvNS_2D():
    det = toyDetector(20000, 0.12, 0.88)
    Na = 6.022e23
    gtoeV = 5.61e32
    mC12 = 12 / Na * gtoeV
    xsec = CEvNS_XS(6, 6, mC12)
    dist = 10
    model0 = SNNuloader(25, "Accretion", "LS220", dist)

    t, _, _, no2d, io2d = calc_visibleSpectra2D(model0, xsec, det, chaname="CEvNS")
    np.savetxt("CEvNS.txt", no2d, fmt="%.3e", delimiter=",")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.imshow(no2d, extent=[-20, 40, 0, 0.1], aspect="auto", cmap="cividis")
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel(r"$E_\nu$ [MeV]", fontsize=14)
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()
    plt.show()


def drawNuE():
    det = toyDetector(20000, 0.12, 0.88)
    xsec = NuE_XS()
    dist = 10
    model0 = SNNuloader(25, "Accretion", "LS220", dist)
    spec = visible_spectrum("NuE", model0, xsec, det)

    t, no, io = calc_visibleSpectra(model0, xsec, det, chaname="NuE")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t*1000, no/1000, "-", lw=2, label=r"$25 M_\odot$, LS220")
    ax.plot(t*1000, io/1000, ":", lw=2, )
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("counts per ms", fontsize=14)
    ax.legend(prop={"size":12})
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()
    plt.savefig(f"timespec_NuE.pdf")
    plt.show()



if __name__ == '__main__':
    drawCEvNS()
    #drawCEvNS_2D()
    #drawNuE()
