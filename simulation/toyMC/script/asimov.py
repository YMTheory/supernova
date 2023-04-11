import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt, cm
from matplotlib import colors
import matplotlib as mpl
import re
plt.style.use("science")
#plt.style.use("dark_background")

def drawTargetDistEthr3D():
    Ethr = np.arange(0, 0.16, 0.01) * 1000
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection='3d')
    colors1 = mpl.cm.viridis(np.linspace(0., 1, 400))

    
    arrIO = np.loadtxt("../results/Garching82703_variedDist_variedTarget_variedEthr_IO_pESeESIBD_MOSensAsimovC14low.csv")
    EthrIO = arrIO[:, 0]
    scaleIO = arrIO[:, 1]
    dT_dataIO_pdfNO = arrIO[:, 2]
    nll_dataIO_pdfNO = arrIO[:, 3]
    dT_dataIO_pdfIO = arrIO[:, 4]
    nll_dataIO_pdfIO = arrIO[:, 5]
    

    for E, sc, nll1, nll2 in zip(EthrIO, scaleIO, nll_dataIO_pdfNO, nll_dataIO_pdfIO):
        idc = int((2*nll1-2*nll2)*10)
        ax.scatter(E*1000, sc*20, 10, s=25, c=colors1[idc], alpha=0.5)


    arrIO = np.loadtxt("../results/Garching82703_variedDist_JUNO_variedEthr_IO_pESeESIBD_MOSensAsimovC14low.csv")
    EthrIO = arrIO[:, 0]
    scaleIO = arrIO[:, 1]
    dT_dataIO_pdfNO = arrIO[:, 2]
    nll_dataIO_pdfNO = arrIO[:, 3]
    dT_dataIO_pdfIO = arrIO[:, 4]
    nll_dataIO_pdfIO = arrIO[:, 5]

    for E, sc, nll1, nll2 in zip(EthrIO, scaleIO, nll_dataIO_pdfNO, nll_dataIO_pdfIO):
        idc = int((2*nll1-2*nll2)*10)
        dist = 10 / np.sqrt(sc)
        ax.scatter(E*1000, 20, dist, s=25, c=colors1[idc], alpha=0.5)


    for Et in range(0, 16, 1):
        Et /= 1000.
        for dist in [5, 6, 7, 8, 9 ,10]:
            for scale in [0.5, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
                NOname = f"../jobs/NOlowC14Ethr{Et:.2f}MeV_target{scale:.1f}_dist{dist}kpc_1D.txt"
                IOname = f"../jobs/IOlowC14Ethr{Et:.2f}MeV_target{scale:.1f}_dist{dist}kpc_1D.txt"
                
                #with open(NOname, "r") as fNO:
                #    for lines in fNO.readlines():
                #        line = lines.strip("\n")
                #        if re.match("Output", line):
                #            data = line.split(" ")
                #            ax.scatter(Et*1000., 20*scale, dist, 2*(float(data[])))
    
                with open(IOname, "r") as fIO:
                    for lines in fIO.readlines():
                        line = lines.strip("\n")
                        if re.match("Output", line):
                            data = line.split(" ")
                            nll1, nll2 = float(data[3]), float(data[5])
                            idc = int(2*(nll1-nll2)*10)
                            print(nll1, nll2, idc)
                            ax.scatter(Et*1000., 20*scale, dist, 2*(nll1-nll2), c=colors1[idc], alpha=0.5)
        




    norm0 = colors.Normalize(0, 40)
    cb1 = plt.colorbar(cm.ScalarMappable(norm=norm0), ax=ax, location="top")
    cb1.set_label(r"$\Delta \chi^2$", fontsize=18)

    ax.set_xlabel(r"$E_\mathrm{thr}$ [keV]", fontsize=16)
    ax.set_zlabel("Distance [kpc]", fontsize=18)
    ax.set_ylabel("Target mass [kton]", fontsize=18)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.tick_params(axis="z", labelsize=16)
    plt.tight_layout()
    plt.show()

    

def drawDeltaTEthr():
    Ethr = np.arange(0, 0.16, 0.01) * 1000

    arrNOnoC14 = np.loadtxt("../jobs/NOnoC14.txt")
    dT_dataNO_pdfNO_noC14_2D = arrNOnoC14[:, 1]
    dT_dataNO_pdfIO_noC14_2D = arrNOnoC14[:, 3]

    arrNOlowC14 = np.loadtxt("../jobs/NOlowC14.txt")
    dT_dataNO_pdfNO_lowC14_2D = arrNOlowC14[:, 1]
    dT_dataNO_pdfIO_lowC14_2D = arrNOlowC14[:, 3]

    arrNOhighC14 = np.loadtxt("../jobs/NOhighC14.txt")
    dT_dataNO_pdfNO_highC14_2D = arrNOhighC14[:, 1]
    dT_dataNO_pdfIO_highC14_2D = arrNOhighC14[:, 3]

    arrIOnoC14 = np.loadtxt("../jobs/IOnoC14.txt")
    dT_dataIO_pdfNO_noC14_2D = arrIOnoC14[:, 1]
    dT_dataIO_pdfIO_noC14_2D = arrIOnoC14[:, 3]

    arrIOlowC14 = np.loadtxt("../jobs/IOlowC14.txt")
    dT_dataIO_pdfNO_lowC14_2D = arrIOlowC14[:, 1]
    dT_dataIO_pdfIO_lowC14_2D = arrIOlowC14[:, 3]

    arrIOhighC14 = np.loadtxt("../jobs/IOhighC14.txt")
    dT_dataIO_pdfNO_highC14_2D = arrIOhighC14[:, 1]
    dT_dataIO_pdfIO_highC14_2D = arrIOhighC14[:, 3]


    arrNOnoC14 = np.loadtxt("../jobs/NOnoC14_ES2D.txt")
    dT_dataNO_pdfNO_noC14_ES2D = arrNOnoC14[:, 1]
    dT_dataNO_pdfIO_noC14_ES2D = arrNOnoC14[:, 3]

    arrNOlowC14 = np.loadtxt("../jobs/NOlowC14_ES2D.txt")
    dT_dataNO_pdfNO_lowC14_ES2D = arrNOlowC14[:, 1]
    dT_dataNO_pdfIO_lowC14_ES2D = arrNOlowC14[:, 3]

    arrNOhighC14 = np.loadtxt("../jobs/NOhighC14_ES2D.txt")
    dT_dataNO_pdfNO_highC14_ES2D = arrNOhighC14[:, 1]
    dT_dataNO_pdfIO_highC14_ES2D = arrNOhighC14[:, 3]

    arrIOnoC14 = np.loadtxt("../jobs/IOnoC14_ES2D.txt")
    dT_dataIO_pdfNO_noC14_ES2D = arrIOnoC14[:, 1]
    dT_dataIO_pdfIO_noC14_ES2D = arrIOnoC14[:, 3]

    arrIOlowC14 = np.loadtxt("../jobs/IOlowC14_ES2D.txt")
    dT_dataIO_pdfNO_lowC14_ES2D = arrIOlowC14[:, 1]
    dT_dataIO_pdfIO_lowC14_ES2D = arrIOlowC14[:, 3]

    arrIOhighC14 = np.loadtxt("../jobs/IOhighC14_ES2D.txt")
    dT_dataIO_pdfNO_highC14_ES2D = arrIOhighC14[:, 1]
    dT_dataIO_pdfIO_highC14_ES2D = arrIOhighC14[:, 3]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(Ethr, dT_dataNO_pdfIO_noC14_2D, "-", color="red", lw=2)
    ax.plot(Ethr, dT_dataNO_pdfIO_noC14_ES2D, ":", color="red", lw=2)
    ax.plot(Ethr, dT_dataNO_pdfIO_lowC14_2D, "-", color="blue", lw=2)
    ax.plot(Ethr, dT_dataNO_pdfIO_lowC14_ES2D, ":", color="blue", lw=2)
    ax.plot(Ethr, dT_dataNO_pdfIO_highC14_2D, "-", color="orange", lw=2)
    ax.plot(Ethr, dT_dataNO_pdfIO_highC14_ES2D, ":", color="orange", lw=2)

    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(prop={"size":10}, loc="upper right", frameon=True, ncol=3)
    ax.set_xlabel(r"$E_\mathrm{thr}$ [keV]", fontsize=18)
    ax.set_ylabel(r"$\Delta t$", fontsize=18)
    ax.tick_params(axis="both", labelsize=16)
    ax.set_title("1D fit", fontsize=18)
    #ax.set_ylim(4, 13)
    
    plt.tight_layout()
    plt.show()



def drawSensEthr():
    Ethr = np.arange(0, 0.16, 0.01) * 1000

    arrNO = np.loadtxt("../results/Garching82703_10kpc_NO_pESeESIBD_MOSensAsimov_JUNO.csv")
    dT_dataNO_pdfNO0 = arrNO[:, 1]
    nll_dataNO_pdfNO0 = arrNO[:, 2]
    dT_dataNO_pdfIO0 = arrNO[:, 3]
    nll_dataNO_pdfIO0 = arrNO[:, 4]
    
    arrIO = np.loadtxt("../results/Garching82703_10kpc_IO_pESeESIBD_MOSensAsimov_JUNO.csv")
    dT_dataIO_pdfNO0 = arrIO[:, 1]
    nll_dataIO_pdfNO0 = arrIO[:, 2]
    dT_dataIO_pdfIO0 = arrIO[:, 3]
    nll_dataIO_pdfIO0 = arrIO[:, 4]
    
    arrNO = np.loadtxt("../results/Garching82703_10kpc_NO_pESeESIBD_MOSensAsimovC14low_JUNO.csv")
    dT_dataNO_pdfNO = arrNO[:, 1]
    nll_dataNO_pdfNO = arrNO[:, 2]
    dT_dataNO_pdfIO = arrNO[:, 3]
    nll_dataNO_pdfIO = arrNO[:, 4]
    
    arrIO = np.loadtxt("../results/Garching82703_10kpc_IO_pESeESIBD_MOSensAsimovC14low_JUNO.csv")
    dT_dataIO_pdfNO = arrIO[:, 1]
    nll_dataIO_pdfNO = arrIO[:, 2]
    dT_dataIO_pdfIO = arrIO[:, 3]
    nll_dataIO_pdfIO = arrIO[:, 4]

    arrNO = np.loadtxt("../results/Garching82703_10kpc_NO_pESeESIBD_MOSensAsimovC14high_JUNO.csv")
    dT_dataNO_pdfNO1 = arrNO[:, 1]
    nll_dataNO_pdfNO1 = arrNO[:, 2]
    dT_dataNO_pdfIO1 = arrNO[:, 3]
    nll_dataNO_pdfIO1 = arrNO[:, 4]
    
    arrIO = np.loadtxt("../results/Garching82703_10kpc_IO_pESeESIBD_MOSensAsimovC14high_JUNO.csv")
    dT_dataIO_pdfNO1 = arrIO[:, 1]
    nll_dataIO_pdfNO1 = arrIO[:, 2]
    dT_dataIO_pdfIO1 = arrIO[:, 3]
    nll_dataIO_pdfIO1 = arrIO[:, 4]
    
    arrNOnoC14 = np.loadtxt("../jobs/NOnoC14.txt")
    nll_dataNO_pdfNO_noC14_2D = arrNOnoC14[:, 2]
    nll_dataNO_pdfIO_noC14_2D = arrNOnoC14[:, 4]

    arrNOlowC14 = np.loadtxt("../jobs/NOlowC14.txt")
    nll_dataNO_pdfNO_lowC14_2D = arrNOlowC14[:, 2]
    nll_dataNO_pdfIO_lowC14_2D = arrNOlowC14[:, 4]

    arrNOhighC14 = np.loadtxt("../jobs/NOhighC14.txt")
    nll_dataNO_pdfNO_highC14_2D = arrNOhighC14[:, 2]
    nll_dataNO_pdfIO_highC14_2D = arrNOhighC14[:, 4]

    arrIOnoC14 = np.loadtxt("../jobs/IOnoC14.txt")
    nll_dataIO_pdfNO_noC14_2D = arrIOnoC14[:, 2]
    nll_dataIO_pdfIO_noC14_2D = arrIOnoC14[:, 4]

    arrIOlowC14 = np.loadtxt("../jobs/IOlowC14.txt")
    nll_dataIO_pdfNO_lowC14_2D = arrIOlowC14[:, 2]
    nll_dataIO_pdfIO_lowC14_2D = arrIOlowC14[:, 4]

    arrIOhighC14 = np.loadtxt("../jobs/IOhighC14.txt")
    nll_dataIO_pdfNO_highC14_2D = arrIOhighC14[:, 2]
    nll_dataIO_pdfIO_highC14_2D = arrIOhighC14[:, 4]


    arrNOnoC14 = np.loadtxt("../jobs/NOnoC14_ES2D.txt")
    nll_dataNO_pdfNO_noC14_ES2D = arrNOnoC14[:, 2]
    nll_dataNO_pdfIO_noC14_ES2D = arrNOnoC14[:, 4]

    arrNOlowC14 = np.loadtxt("../jobs/NOlowC14_ES2D.txt")
    nll_dataNO_pdfNO_lowC14_ES2D = arrNOlowC14[:, 2]
    nll_dataNO_pdfIO_lowC14_ES2D = arrNOlowC14[:, 4]

    arrNOhighC14 = np.loadtxt("../jobs/NOhighC14_ES2D.txt")
    nll_dataNO_pdfNO_highC14_ES2D = arrNOhighC14[:, 2]
    nll_dataNO_pdfIO_highC14_ES2D = arrNOhighC14[:, 4]

    arrIOnoC14 = np.loadtxt("../jobs/IOnoC14_ES2D.txt")
    nll_dataIO_pdfNO_noC14_ES2D = arrIOnoC14[:, 2]
    nll_dataIO_pdfIO_noC14_ES2D = arrIOnoC14[:, 4]

    arrIOlowC14 = np.loadtxt("../jobs/IOlowC14_ES2D.txt")
    nll_dataIO_pdfNO_lowC14_ES2D = arrIOlowC14[:, 2]
    nll_dataIO_pdfIO_lowC14_ES2D = arrIOlowC14[:, 4]

    arrIOhighC14 = np.loadtxt("../jobs/IOhighC14_ES2D.txt")
    nll_dataIO_pdfNO_highC14_ES2D = arrIOhighC14[:, 2]
    nll_dataIO_pdfIO_highC14_ES2D = arrIOhighC14[:, 4]





    """
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(Ethr, 2*(nll_dataIO_pdfNO0 - nll_dataIO_pdfIO0), "-", lw=3, color="crimson",    label=r"IO, 1D fit")
    ax.plot(Ethr, 2*(nll_dataNO_pdfIO0 - nll_dataNO_pdfNO0), "-", lw=3, color="blue",       label=r"NO, 1D fit")
    ax.plot(Ethr, 2*(nll_dataIO_pdfNO_noC14_2D - nll_dataIO_pdfIO_noC14_2D), "--", lw=3, color="crimson", label=r"IO, 2D fit")
    ax.plot(Ethr, 2*(nll_dataNO_pdfIO_noC14_2D - nll_dataNO_pdfNO_noC14_2D), "--", lw=3, color="blue", label=r"NO, 2D fit")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(prop={"size":12}, loc="upper right", frameon=True, ncol=2)
    ax.set_xlabel(r"$E_\mathrm{thr}$ [keV]", fontsize=18)
    ax.set_ylabel(r"$\Delta \chi^2$", fontsize=18)
    ax.tick_params(axis="both", labelsize=16)
    ax.set_ylim(6, 13)
    plt.tight_layout()
    plt.savefig("../plots/Garching82703_10kpc_JUNO_AsimovMOSensEthr_fit1Dtofit2D.pdf")
    plt.show()
    
    """
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    
    ax[0].plot(Ethr, 2*(nll_dataIO_pdfNO0 - nll_dataIO_pdfIO0), "-", lw=3, color="crimson",    label=r"IO w/o $^{14}$C")
    ax[0].plot(Ethr, 2*(nll_dataNO_pdfIO0 - nll_dataNO_pdfNO0), "-", lw=3, color="blue",       label=r"NO w/o $^{14}$C")
    ax[0].plot(Ethr, 2*(nll_dataIO_pdfNO - nll_dataIO_pdfIO),  "--", lw=3, color="crimson",    label=r"$1\times10^{-18}$ g/g")
    ax[0].plot(Ethr, 2*(nll_dataNO_pdfIO - nll_dataNO_pdfNO),  "--", lw=3, color="blue",       label=r"$1\times10^{-18}$ g/g")
    ax[0].plot(Ethr, 2*(nll_dataIO_pdfNO1 - nll_dataIO_pdfIO1), ":", lw=3, color="crimson",    label=r"$1\times10^{-17}$ g/g")
    ax[0].plot(Ethr, 2*(nll_dataNO_pdfIO1 - nll_dataNO_pdfNO1), ":", lw=3, color="blue",       label=r"$1\times10^{-17}$ g/g")
    #ax.fill_between(Ethr, 2*nll_dataIO_pdfNO - 2*nll_dataIO_pdfIO, 2*nll_dataIO_pdfNO1 - 2*nll_dataIO_pdfIO1, color="crimson", hatch="x", alpha=0.5, label="IO, 1D-fit")
            #label=r"IO w/ $^{14}$C from $1\times10^{-17}$ g/g -$1\times10^{-18}$ g/g")
    #ax.fill_between(Ethr, 2*nll_dataNO_pdfIO - 2*nll_dataNO_pdfNO, 2*nll_dataNO_pdfIO1 - 2*nll_dataNO_pdfNO1, color="blue",    hatch="x", alpha=0.5, label="NO, 1D-fit")
            #label=r"NO w/ $^{14}$C from $1\times10^{-17}$ g/g -$1\times10^{-18}$ g/g")
    # NO: noC14 2D fit: Output 82703 0.02 -0.9712528682742865 74.44348998451241 -0.0020660594210630385 78.6544439093938
    # IO: noC14 2D fit: Output 82703 0.02 1.1923428113065802 79.24338286435079 -0.0009695025772358991 74.38532596510849
    # ax.plot(100, 2*(4.094906049170585 + 0.04801026296706014), "o", ms=8, fillstyle="none", color="blue")
    # ax.plot(100, 2*(4.597644204572177 + 0.10503156100216483), "s", ms=8, fillstyle="none", color="crimson")
    # ax.plot(100, 2*(78.6544439093938 - 74.44348998451241), "p", ms=8, fillstyle="none", color="blue")
    # ax.plot(100, 2*(79.24338286435079 - 74.38532596510849), "d", ms=8, fillstyle="none", color="crimson")
    ax[1].plot(Ethr, 2*(nll_dataIO_pdfNO_noC14_2D - nll_dataIO_pdfIO_noC14_2D), "-", lw=3, color="orange", label=r"IO w/o $^{14}$C")
    ax[1].plot(Ethr, 2*(nll_dataNO_pdfIO_noC14_2D - nll_dataNO_pdfNO_noC14_2D), "-", lw=3, color="green", label=r"NO w/o $^{14}$C")
    ax[1].plot(Ethr, 2*(nll_dataIO_pdfNO_lowC14_2D - nll_dataIO_pdfIO_lowC14_2D),  "--", lw=3, color="orange", label=r"$1\times10^{-18}$ g/g")
    ax[1].plot(Ethr, 2*(nll_dataNO_pdfIO_lowC14_2D - nll_dataNO_pdfNO_lowC14_2D),  "--", lw=3, color="green", label=r"$1\times10^{-18}$ g/g")
    ax[1].plot(Ethr, 2*(nll_dataIO_pdfNO_highC14_2D - nll_dataIO_pdfIO_highC14_2D), ":", lw=3, color="orange", label=r"$1\times10^{-17}$ g/g")
    ax[1].plot(Ethr, 2*(nll_dataNO_pdfIO_highC14_2D - nll_dataNO_pdfNO_highC14_2D), ":", lw=3, color="green", label=r"$1\times10^{-17}$ g/g")


    ax[2].plot(Ethr, 2*(nll_dataIO_pdfNO_noC14_ES2D - nll_dataIO_pdfIO_noC14_ES2D), "-", lw=3, color="orange", label=r"IO w/o $^{14}$C")
    ax[2].plot(Ethr, 2*(nll_dataNO_pdfIO_noC14_ES2D - nll_dataNO_pdfNO_noC14_ES2D), "-", lw=3, color="green", label=r"NO w/o $^{14}$C")
    ax[2].plot(Ethr, 2*(nll_dataIO_pdfNO_lowC14_ES2D - nll_dataIO_pdfIO_lowC14_ES2D),  "--", lw=3, color="orange", label=r"$1\times10^{-18}$ g/g")
    ax[2].plot(Ethr, 2*(nll_dataNO_pdfIO_lowC14_ES2D - nll_dataNO_pdfNO_lowC14_ES2D),  "--", lw=3, color="green", label=r"$1\times10^{-18}$ g/g")
    ax[2].plot(Ethr, 2*(nll_dataIO_pdfNO_highC14_ES2D - nll_dataIO_pdfIO_highC14_ES2D), ":", lw=3, color="orange", label=r"$1\times10^{-17}$ g/g")
    ax[2].plot(Ethr, 2*(nll_dataNO_pdfIO_highC14_ES2D - nll_dataNO_pdfNO_highC14_ES2D), ":", lw=3, color="green", label=r"$1\times10^{-17}$ g/g")


    #ax.fill_between(Ethr, 2*nll_dataIO_pdfNO_lowC14_2D - 2*nll_dataIO_pdfIO_lowC14_2D, 2*nll_dataIO_pdfNO_highC14_2D - 2*nll_dataIO_pdfIO_highC14_2D, hatch="|", color="crimson", alpha=0.5, label=r"IO, 2D fit")
    #ax.fill_between(Ethr, 2*nll_dataNO_pdfIO_lowC14_2D - 2*nll_dataNO_pdfNO_lowC14_2D, 2*nll_dataNO_pdfIO_highC14_2D - 2*nll_dataNO_pdfNO_highC14_2D, hatch="|", color="blue", alpha=0.5, label=r"NO, 2D fit")

    #ax.plot(Ethr, 2 * (nll_dataNO_pdfIO_noC14_2D - nll_dataNO_pdfNO_noC14_2D),   ":", lw=3, label="NO, 2D, noC14")
    #ax.plot(Ethr, -2 * (nll_dataIO_pdfIO_noC14_2D - nll_dataIO_pdfNO_noC14_2D),   ":", lw=2, label="IO, 2D, noC14")
    #ax.plot(Ethr, 2 * (nll_dataNO_pdfIO_lowC14_2D - nll_dataNO_pdfNO_lowC14_2D), ":",  color="blue",    lw=3, label="NO, 2D, lowC14")
    #ax.plot(Ethr, -2 * (nll_dataIO_pdfIO_lowC14_2D - nll_dataIO_pdfNO_lowC14_2D), ":", color="crimson", lw=3, label="IO, 2D, lowC14")

    ax[0].grid(True, linestyle=":", alpha=0.5)
    #ax[0].legend(prop={"size":10}, loc="upper right", frameon=True, ncol=3)
    ax[0].set_xlabel(r"$E_\mathrm{thr}$ [keV]", fontsize=18)
    ax[0].set_ylabel(r"$\Delta \chi^2$", fontsize=18)
    ax[0].tick_params(axis="both", labelsize=16)
    ax[0].set_title("1D fit", fontsize=18)
    ax[0].set_ylim(4, 13)
    ax[1].grid(True, linestyle=":", alpha=0.5)
    #ax[1].legend(prop={"size":10}, loc="lower center",frameon=True, ncol=3)
    ax[1].set_xlabel(r"$E_\mathrm{thr}$ [keV]", fontsize=18)
    ax[1].set_ylabel(r"$\Delta \chi^2$", fontsize=18)
    ax[1].tick_params(axis="both", labelsize=16)
    ax[1].set_title("pES 2D fit", fontsize=18)
    ax[1].set_ylim(4, 13)
    ax[2].grid(True, linestyle=":", alpha=0.5)
    #ax[2].legend(prop={"size":10}, loc="lower center",frameon=True, ncol=3)
    ax[2].set_xlabel(r"$E_\mathrm{thr}$ [keV]", fontsize=18)
    ax[2].set_ylabel(r"$\Delta \chi^2$", fontsize=18)
    ax[2].tick_params(axis="both", labelsize=16)
    ax[2].set_title("pES+eES 2D fit", fontsize=18)
    ax[2].set_ylim(4, 13)

    plt.subplots_adjust(hspace=0.0, wspace=0.2, bottom=0.16, left=0.05, right=0.95)
    #plt.savefig("../plots/AsimovTest_EthrSensC14low2high_JUNO_.pdf")
    plt.show()

def drawSensTimeEthr():
    Ethr = np.arange(0, 0.16, 0.01)

    arrNO = np.loadtxt("../results/Garching82703_10kpc_NO_pESeESIBD_MOSensAsimov_JUNO.csv")
    dT_dataNO_pdfNO0 = arrNO[:, 1]
    nll_dataNO_pdfNO0 = arrNO[:, 2]
    dT_dataNO_pdfIO0 = arrNO[:, 3]
    nll_dataNO_pdfIO0 = arrNO[:, 4]
    
    arrIO = np.loadtxt("../results/Garching82703_10kpc_IO_pESeESIBD_MOSensAsimov_JUNO.csv")
    dT_dataIO_pdfNO0 = arrIO[:, 1]
    nll_dataIO_pdfNO0 = arrIO[:, 2]
    dT_dataIO_pdfIO0 = arrIO[:, 3]
    nll_dataIO_pdfIO0 = arrIO[:, 4]
    
    arrNO = np.loadtxt("../results/Garching82703_10kpc_NO_pESeESIBD_MOSensAsimovC14low_JUNO.csv")
    dT_dataNO_pdfNO = arrNO[:, 1]
    nll_dataNO_pdfNO = arrNO[:, 2]
    dT_dataNO_pdfIO = arrNO[:, 3]
    nll_dataNO_pdfIO = arrNO[:, 4]
    
    arrIO = np.loadtxt("../results/Garching82703_10kpc_IO_pESeESIBD_MOSensAsimovC14low_JUNO.csv")
    dT_dataIO_pdfNO = arrIO[:, 1]
    nll_dataIO_pdfNO = arrIO[:, 2]
    dT_dataIO_pdfIO = arrIO[:, 3]
    nll_dataIO_pdfIO = arrIO[:, 4]

    arrNO = np.loadtxt("../results/Garching82703_10kpc_NO_pESeESIBD_MOSensAsimovC14high_JUNO.csv")
    dT_dataNO_pdfNO1 = arrNO[:, 1]
    nll_dataNO_pdfNO1 = arrNO[:, 2]
    dT_dataNO_pdfIO1 = arrNO[:, 3]
    nll_dataNO_pdfIO1 = arrNO[:, 4]
    
    arrIO = np.loadtxt("../results/Garching82703_10kpc_IO_pESeESIBD_MOSensAsimovC14high_JUNO.csv")
    dT_dataIO_pdfNO1 = arrIO[:, 1]
    nll_dataIO_pdfNO1 = arrIO[:, 2]
    dT_dataIO_pdfIO1 = arrIO[:, 3]
    nll_dataIO_pdfIO1 = arrIO[:, 4]
    
    
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    
    ax[0].plot(Ethr, dT_dataNO_pdfIO0 - dT_dataNO_pdfNO0, "-", lw=2, color="cyan")
    ax[0].fill_between(Ethr, dT_dataNO_pdfIO - dT_dataNO_pdfNO, dT_dataNO_pdfIO1 - dT_dataNO_pdfNO1, color="cyan", alpha=0.5)
    ax[0].set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax[0].set_ylabel(r"$\Delta t(\mathrm{IO}) - \Delta t(\mathrm{NO})$ [ms]", fontsize=16, color="cyan")
    ax[0].tick_params(axis="y", labelsize=15, labelcolor="cyan")
    ax[0].tick_params(axis="x", labelsize=15, labelcolor="black")
    ax[0].set_ylim(-5, 0)
    ax[0].set_title("NO data", fontsize=16)
    ax[0].grid(True, linestyle=":", alpha=0.5)
    
    ax0 = ax[0].twinx()
    ax0.plot(Ethr, 2*(nll_dataNO_pdfIO0 - nll_dataNO_pdfNO0), "-", lw=3, color="violet")
    ax0.fill_between(Ethr, 2*nll_dataNO_pdfIO - 2*nll_dataNO_pdfNO, 2*nll_dataNO_pdfIO1 - 2*nll_dataNO_pdfNO1, color="violet", alpha=0.5)
    ax0.fill_between(Ethr, 100+2*nll_dataNO_pdfIO - 2*nll_dataNO_pdfNO, 100+2*nll_dataNO_pdfIO1 - 2*nll_dataNO_pdfNO1, color="gray", alpha=0.5, label=r"$^{14}$C: $1\times10{-18}$ g/g - $1\times10^{-17}$ g/g")
    ax0.plot(Ethr, 100+2*(nll_dataNO_pdfIO0 - nll_dataNO_pdfNO0), "-", lw=3, color="gray", label=r"no $^{14}$C")
    ax0.set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax0.set_ylabel(r"$\Delta \chi^2$", fontsize=16, color="violet", rotation=270, labelpad=20)
    ax0.tick_params(axis="y", labelsize=15, labelcolor="violet")
    ax0.set_ylim(3, 10)
    ax0.legend(prop={"size":12}, loc="lower left", frameon=True)
    
    
    ax[1].plot(Ethr, dT_dataIO_pdfIO0 - dT_dataIO_pdfNO0, "-", lw=2, color="cyan", label=r"$\Delta t$")
    ax[1].plot(Ethr, 200+dT_dataIO_pdfIO0 - dT_dataIO_pdfNO0, "-", lw=2, color="violet", label=r"$\Delta \chi^2$")
    ax[1].set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax[1].fill_between(Ethr, dT_dataIO_pdfIO - dT_dataIO_pdfNO, dT_dataIO_pdfIO1 - dT_dataIO_pdfNO1, color="cyan", alpha=0.5)
    ax[1].set_ylabel(r"$\Delta t(\mathrm{IO}) - \Delta t(\mathrm{NO})$ [ms]", fontsize=16, color="cyan")
    ax[1].tick_params(axis="y", labelsize=15, labelcolor="cyan")
    ax[1].tick_params(axis="x", labelsize=15, labelcolor="black")
    ax[1].set_ylim(-5, 0)
    ax[1].set_title("IO data", fontsize=16)
    ax[1].legend(prop={"size":12}, loc="lower right", frameon=False)
    
    ax1 = ax[1].twinx()
    ax1.plot(Ethr, 2*(nll_dataIO_pdfNO0 - nll_dataIO_pdfIO0), "-", lw=3, color="violet")
    ax1.fill_between(Ethr, 2*nll_dataIO_pdfNO - 2*nll_dataIO_pdfIO, 2*nll_dataIO_pdfNO1 - 2*nll_dataIO_pdfIO1, color="violet", alpha=0.5)
    ax1.set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax1.set_ylabel(r"$\Delta \chi^2$", fontsize=16, color="violet", rotation=270, labelpad=20)
    ax1.tick_params(axis="y", labelsize=15, labelcolor="violet")
    ax1.set_ylim(3, 12)
    ax[1].grid(True, linestyle=":", alpha=0.5)
    
    plt.subplots_adjust(hspace=0.0, wspace=0.4, bottom=0.16)
    plt.savefig("../plots/AsimovTest_EthrTimeSensC14low2high_JUNO_dark.pdf")
    plt.show()


def drawSensTargetEthr():

    arrNO = np.loadtxt("../results/Garching82703_variedDist_variedTarget_variedEthr_NO_pESeESIBD_MOSensAsimovC14low.csv")
    EthrNO = arrNO[:, 0]
    scaleNO = arrNO[:, 1]
    dT_dataNO_pdfNO = arrNO[:, 2]
    nll_dataNO_pdfNO = arrNO[:, 3]
    dT_dataNO_pdfIO = arrNO[:, 4]
    nll_dataNO_pdfIO = arrNO[:, 5]
    
    arrIO = np.loadtxt("../results/Garching82703_variedDist_variedTarget_variedEthr_IO_pESeESIBD_MOSensAsimovC14low.csv")
    EthrIO = arrIO[:, 0]
    scaleIO = arrIO[:, 1]
    dT_dataIO_pdfNO = arrIO[:, 2]
    nll_dataIO_pdfNO = arrIO[:, 3]
    dT_dataIO_pdfIO = arrIO[:, 4]
    nll_dataIO_pdfIO = arrIO[:, 5]
    
    imNO = np.zeros((6, 8))  # x-axis: scale; y-axis: Ethr
    imIO = np.zeros((6, 8))
    for E, s, nll1, nll2 in zip(EthrNO, scaleNO, nll_dataNO_pdfNO, nll_dataNO_pdfIO):
        idy = int(s/0.5) - 1
        idx = 5 - round((E-0.10)*100)
        imNO[idx, idy] = 2 * (nll2 - nll1)
    for E, s, nll1, nll2 in zip(EthrIO, scaleIO, nll_dataIO_pdfNO, nll_dataIO_pdfIO):
        idy = int(s/0.5) - 1
        idx = 5 - round((E-0.10)*100)
        imIO[idx, idy] = -2 * (nll2 - nll1)


    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    im0 = ax[0].imshow(imNO, extent=[5, 85, 95, 155], aspect="auto")
    ax[0].set_xlabel(r"Target mass [kton]", fontsize=16)
    ax[0].set_ylabel(r"$E_\mathrm{thr}$ [keV]", fontsize=16)
    ax[0].tick_params(axis="y", labelsize=15)
    ax[0].tick_params(axis="x", labelsize=15)
    ax[0].set_title("NO data", fontsize=16)
    cb0 = plt.colorbar(im0, ax=ax[0])
    cb0.set_label(r"$\Delta \chi^2$", fontsize=16, rotation=270, labelpad=25)
    cb0.ax.tick_params(labelsize=15)
    X0, Y0 = np.arange(10, 90, 10), np.arange(100, 160, 10)
    cs0 = ax[0].contour(X0, Y0, imNO, levels=[9, 16, 25], colors=["red", "darkorange", "darkviolet"], linewidth=2)
    ax[0].text(0.6*20, 150, r"$\Delta\chi^2=9$", fontsize=15, color="red")
    ax[0].text(1.9*20, 150, r"$\Delta\chi^2=16$", fontsize=15, color="darkorange")
    ax[0].text(3.2*20, 150, r"$\Delta\chi^2=25$", fontsize=15, color="darkviolet")
    ax[0].plot(1*20, 100., "*", ms=10, color="chocolate")
    ax[0].text(1.2*20, 98, r"$\Delta \chi^2 = 6.7$ (JUNO)", color="chocolate", fontsize=15)
    ax[0].set_xticks(np.arange(10, 90, 10), ["10", "20", "30", "40", "50", "60", "70", "80"])

    
    im1 = ax[1].imshow(imIO, extent=[5, 85, 95, 155], aspect="auto")
    ax[1].set_xlabel(r"Target mass [kton]", fontsize=16)
    ax[1].set_ylabel(r"$E_\mathrm{thr}$ [keV]", fontsize=16)
    ax[1].tick_params(axis="y", labelsize=15)
    ax[1].tick_params(axis="x", labelsize=15)
    ax[1].set_title("IO data", fontsize=16)
    cb1 = plt.colorbar(im1, ax=ax[1])
    cb1.set_label(r"$\Delta \chi^2$", fontsize=16, rotation=270, labelpad=25)
    cb1.ax.tick_params(labelsize=15)
    X1, Y1 = np.arange(10, 90, 10), np.arange(100, 160, 10)
    cs1 = ax[1].contour(X1, Y1, imIO, levels=[9, 16, 25], colors=["red", "darkorange", "darkviolet"], linewidth=2)
    ax[1].text(0.6*20, 150., r"$\Delta\chi^2=9$", fontsize=15, color="red")
    ax[1].text(1.9*20, 150., r"$\Delta\chi^2=16$", fontsize=15, color="darkorange")
    ax[1].text(3.2*20, 150., r"$\Delta\chi^2=25$", fontsize=15, color="darkviolet")
    ax[1].plot(1*20, 100., "*", ms=10, color="chocolate")
    ax[1].text(1.2*20, 98, r"$\Delta \chi^2 = 7.4$ (JUNO)", color="chocolate", fontsize=15)
    ax[1].set_xticks(np.arange(10, 90, 10), ["10", "20", "30", "40", "50", "60", "70", "80"])
    
    plt.tight_layout()
    plt.savefig("../plots/AsimovTest_TargetEthrSensC14low_dark.pdf")
    plt.show()


def drawSensDistEthr():

    arrNO = np.loadtxt("../results/Garching82703_variedDist_JUNO_variedEthr_NO_pESeESIBD_MOSensAsimovC14low.csv")
    EthrNO = arrNO[:, 0]
    scaleNO = arrNO[:, 1]
    dT_dataNO_pdfNO = arrNO[:, 2]
    nll_dataNO_pdfNO = arrNO[:, 3]
    dT_dataNO_pdfIO = arrNO[:, 4]
    nll_dataNO_pdfIO = arrNO[:, 5]
    
    arrIO = np.loadtxt("../results/Garching82703_variedDist_JUNO_variedEthr_IO_pESeESIBD_MOSensAsimovC14low.csv")
    EthrIO = arrIO[:, 0]
    scaleIO = arrIO[:, 1]
    dT_dataIO_pdfNO = arrIO[:, 2]
    nll_dataIO_pdfNO = arrIO[:, 3]
    dT_dataIO_pdfIO = arrIO[:, 4]
    nll_dataIO_pdfIO = arrIO[:, 5]
    
    imNO = np.zeros((6, 6))  # x-axis: scale; y-axis: Ethr
    imIO = np.zeros((6, 6))
    for E, s, nll1, nll2 in zip(EthrNO, scaleNO, nll_dataNO_pdfNO, nll_dataNO_pdfIO):
        d = int(np.sqrt(100/s))
        idy = d - 5
        idx = 5 - round((E-0.10)*100)
        imNO[idx, idy] = 2 * (nll2 - nll1)
    for E, s, nll1, nll2 in zip(EthrIO, scaleIO, nll_dataIO_pdfNO, nll_dataIO_pdfIO):
        d = int(np.sqrt(100/s))
        idy = d - 5
        idx = 5 - round((E-0.10)*100)
        imIO[idx, idy] = -2 * (nll2 - nll1)


    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    im0 = ax[0].imshow(imNO, extent=[4.5, 10.5, 95, 155], aspect="auto")
    ax[0].set_xlabel(r"Distance [kpc]", fontsize=16)
    ax[0].set_ylabel(r"$E_\mathrm{thr}$ [keV]", fontsize=16)
    ax[0].tick_params(axis="y", labelsize=15)
    ax[0].tick_params(axis="x", labelsize=15)
    ax[0].set_title("NO data", fontsize=16)
    cb0 = plt.colorbar(im0, ax=ax[0])
    cb0.set_label(r"$\Delta \chi^2$", fontsize=16, rotation=270, labelpad=25)
    cb0.ax.tick_params(labelsize=15)
    X0, Y0 = np.arange(5, 11, 1), np.arange(100, 160, 10)
    cs0 = ax[0].contour(X0, Y0, imNO, levels=[9, 16, 25], colors=["red", "darkorange", "darkviolet"], linewidth=2)
    ax[0].text(8.7, 150, r"$\Delta\chi^2=9$", fontsize=15, color="red")
    ax[0].text(6.9, 150, r"$\Delta\chi^2=16$", fontsize=15, color="darkorange")
    ax[0].text(5, 150, r"$\Delta\chi^2=25$", fontsize=15, color="darkviolet")
    ax[0].set_xticks(np.arange(5, 11, 1), ["5", "6", "7", "8", "9", "10"])

    
    im1 = ax[1].imshow(imIO, extent=[4.5, 10.5, 95, 155], aspect="auto")
    ax[1].set_xlabel(r"Distance [kpc]", fontsize=16)
    ax[1].set_ylabel(r"$E_\mathrm{thr}$ [keV]", fontsize=16)
    ax[1].tick_params(axis="y", labelsize=15)
    ax[1].tick_params(axis="x", labelsize=15)
    ax[1].set_title("IO data", fontsize=16)
    cb1 = plt.colorbar(im1, ax=ax[1])
    cb1.set_label(r"$\Delta \chi^2$", fontsize=16, rotation=270, labelpad=25)
    cb1.ax.tick_params(labelsize=15)
    X1, Y1 = np.arange(5, 11, 1), np.arange(100, 160, 10)
    cs0 = ax[1].contour(X1, Y1, imIO, levels=[9, 16, 25], colors=["red", "darkorange", "darkviolet"], linewidth=2)
    ax[1].text(8.7, 150, r"$\Delta\chi^2=9$", fontsize=15, color="red")
    ax[1].text(6.9, 150, r"$\Delta\chi^2=16$", fontsize=15, color="darkorange")
    ax[1].text(5, 150, r"$\Delta\chi^2=25$", fontsize=15, color="darkviolet")
    ax[1].set_xticks(np.arange(5, 11, 1), ["5", "6", "7", "8", "9", "10"])
    
    plt.tight_layout()
    plt.savefig("../plots/AsimovTest_DistEthrSensC14low_dark.pdf")
    plt.show()

def compareExps():
    
    dist = [10, 9, 8, 7, 6, 5]

    arrNO_JUNO = np.loadtxt("../results/Garching82703_variedDist_JUNO_NO_pESeESIBD_MOSensAsimov.csv")
    EthrNO_JUNO = arrNO_JUNO[:, 0]
    nll_dataNO_pdfNO_JUNO = arrNO_JUNO[:, 3]
    nll_dataNO_pdfIO_JUNO = arrNO_JUNO[:, 5]

    arrIO_JUNO = np.loadtxt("../results/Garching82703_variedDist_JUNO_IO_pESeESIBD_MOSensAsimov.csv")
    EthrIO_JUNO = arrIO_JUNO[:, 0]
    nll_dataIO_pdfNO_JUNO = arrIO_JUNO[:, 3]
    nll_dataIO_pdfIO_JUNO = arrIO_JUNO[:, 5]

    arrNO_JUNO0 = np.loadtxt("../results/Garching82703_variedDist_JUNO_NO_pESeESshort_MOSensAsimov.csv")
    EthrNO_JUNO0 = arrNO_JUNO0[:, 0]
    nll_dataNO_pdfNO_JUNO0 = arrNO_JUNO0[:, 3]
    nll_dataNO_pdfIO_JUNO0 = arrNO_JUNO0[:, 5]

    arrIO_JUNO0 = np.loadtxt("../results/Garching82703_variedDist_JUNO_IO_pESeESshort_MOSensAsimov.csv")
    EthrIO_JUNO0 = arrIO_JUNO0[:, 0]
    nll_dataIO_pdfNO_JUNO0 = arrIO_JUNO0[:, 3]
    nll_dataIO_pdfIO_JUNO0 = arrIO_JUNO0[:, 5]

    arrNO_THEIA = np.loadtxt("../results/Garching82703_variedDist_THEIA100_NO_pESeES_MOSensAsimov.csv")
    EthrNO_THEIA = arrNO_THEIA[:, 0]
    nll_dataNO_pdfNO_THEIA = arrNO_THEIA[:, 3]
    nll_dataNO_pdfIO_THEIA = arrNO_THEIA[:, 5]

    arrIO_THEIA = np.loadtxt("../results/Garching82703_variedDist_THEIA100_IO_pESeES_MOSensAsimov.csv")
    EthrIO_THEIA = arrIO_THEIA[:, 0]
    nll_dataIO_pdfNO_THEIA = arrIO_THEIA[:, 3]
    nll_dataIO_pdfIO_THEIA = arrIO_THEIA[:, 5]

    arrNO_THEIA25 = np.loadtxt("../results/Garching82703_variedDist_THEIA25_NO_pESeESshort_MOSensAsimov.csv")
    EthrNO_THEIA25 = arrNO_THEIA25[:, 0]
    nll_dataNO_pdfNO_THEIA25 = arrNO_THEIA25[:, 3]
    nll_dataNO_pdfIO_THEIA25 = arrNO_THEIA25[:, 5]

    arrIO_THEIA25 = np.loadtxt("../results/Garching82703_variedDist_THEIA25_IO_pESeESshort_MOSensAsimov.csv")
    EthrIO_THEIA25 = arrIO_THEIA25[:, 0]
    nll_dataIO_pdfNO_THEIA25 = arrIO_THEIA25[:, 3]
    nll_dataIO_pdfIO_THEIA25 = arrIO_THEIA25[:, 5]
    
    arrNO_THEIA25all = np.loadtxt("../results/Garching82703_variedDist_THEIA25_NO_pESeESIBD_MOSensAsimov.csv")
    EthrNO_THEIA25all = arrNO_THEIA25all[:, 0]
    nll_dataNO_pdfNO_THEIA25all = arrNO_THEIA25all[:, 3]
    nll_dataNO_pdfIO_THEIA25all = arrNO_THEIA25all[:, 5]

    arrIO_THEIA25all = np.loadtxt("../results/Garching82703_variedDist_THEIA25_IO_pESeESIBD_MOSensAsimov.csv")
    EthrIO_THEIA25all = arrIO_THEIA25all[:, 0]
    nll_dataIO_pdfNO_THEIA25all = arrIO_THEIA25all[:, 3]
    nll_dataIO_pdfIO_THEIA25all = arrIO_THEIA25all[:, 5]
    
    
    arrNO_HK = np.loadtxt("../results/Garching82703_variedDist_HyperK_NO_eESshort_MOSensAsimov.csv")
    EthrNO_HK = arrNO_HK[:, 0]
    nll_dataNO_pdfNO_HK = arrNO_HK[:, 3]
    nll_dataNO_pdfIO_HK = arrNO_HK[:, 5]

    arrIO_HK = np.loadtxt("../results/Garching82703_variedDist_HyperK_IO_eESshort_MOSensAsimov.csv")
    EthrIO_HK = arrIO_HK[:, 0]
    nll_dataIO_pdfNO_HK = arrIO_HK[:, 3]
    nll_dataIO_pdfIO_HK = arrIO_HK[:, 5]
   
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(dist, 2 * (nll_dataNO_pdfIO_JUNO - nll_dataNO_pdfNO_JUNO),   "-", color="royalblue",  lw=2.5, label="JUNO: default")
    ax[0].plot(dist, 2 * (nll_dataNO_pdfIO_JUNO0 - nll_dataNO_pdfNO_JUNO0), "--",color="royalblue",  lw=2.5, label="JUNO: ES-only")
    ax[0].plot(dist, 2 * (nll_dataNO_pdfIO_THEIA25all - nll_dataNO_pdfNO_THEIA25all), "-",color="violet",  lw=2.5, label="THEIA-25: default")
    ax[0].plot(dist, 2 * (nll_dataNO_pdfIO_THEIA25 - nll_dataNO_pdfNO_THEIA25), "--",color="violet",  lw=2.5, label="THEIA-25: ES-only")
    ax[0].plot(dist, 2 * (nll_dataNO_pdfIO_THEIA - nll_dataNO_pdfNO_THEIA), "-.",color="orange",  lw=2.5, label="THEIA-100: ES-only")
    ax[0].plot(dist, 2 * (nll_dataNO_pdfIO_HK - nll_dataNO_pdfNO_HK),  ":",lw=2.5, color="green", label="HyperK: ES-only")
    ax[0].set_xlabel(r"Distance [kpc]", fontsize=16)
    ax[0].set_ylabel(r"$\Delta \chi^2$", fontsize=16)
    ax[0].tick_params(axis="y", labelsize=15)
    ax[0].tick_params(axis="x", labelsize=15)
    ax[0].set_title("NO data", fontsize=16)
    #ax[0].legend(prop={"size":13})
    ax[0].set_ylim(0, 40)
    ax[0].grid(True, linestyle=":")

    ax[1].plot(dist, -2 * (nll_dataIO_pdfIO_JUNO - nll_dataIO_pdfNO_JUNO),"-", color="royalblue", lw=2.5, label="JUNO: default")
    ax[1].plot(dist, -2 * (nll_dataIO_pdfIO_JUNO0 - nll_dataIO_pdfNO_JUNO0),"--", color="royalblue", lw=2.5, label="JUNO: ES only")
    ax[1].plot(dist, -2 * (nll_dataIO_pdfIO_THEIA25all - nll_dataIO_pdfNO_THEIA25all), "-", color="violet", lw=2.5, label="THEIA-25: default")
    ax[1].plot(dist, -2 * (nll_dataIO_pdfIO_THEIA25 - nll_dataIO_pdfNO_THEIA25), "--", color="violet", lw=2.5, label="THEIA-25: ES only")
    ax[1].plot(dist, -2 * (nll_dataIO_pdfIO_THEIA - nll_dataIO_pdfNO_THEIA), "-.", color="orange", lw=2.5, label="THEIA-100: ES only")
    ax[1].plot(dist,  -2 * (nll_dataIO_pdfIO_HK - nll_dataIO_pdfNO_HK), ":",lw=2.5, color="green", label="HyperK: ES only")
    ax[1].set_xlabel(r"Distance [kpc]", fontsize=16)
    ax[1].set_ylabel(r"$\Delta \chi^2$", fontsize=16)
    ax[1].tick_params(axis="y", labelsize=15)
    ax[1].tick_params(axis="x", labelsize=15)
    ax[1].legend(prop={"size":13}, frameon=True)
    ax[1].set_title("IO data", fontsize=16)
    ax[1].set_ylim(0, 40)
    ax[1].grid(True, linestyle=":")

    plt.tight_layout()
    plt.savefig("../plots/compareExps_dark.pdf")
    plt.show()



def compareModels():

    arrNO_Bur = np.loadtxt("../results/Burrows_10kpc_NO_noC14_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modNO_Bur = arrNO_Bur[:, 0]
    tmaxNO_Bur = arrNO_Bur[:, 1]
    nll_dataNO_pdfNO_Bur = arrNO_Bur[:, 3]
    nll_dataNO_pdfIO_Bur = arrNO_Bur[:, 5]

    arrIO_Bur = np.loadtxt("../results/Burrows_10kpc_IO_noC14_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modIO_Bur = arrIO_Bur[:, 0]
    tmaxIO_Bur = arrIO_Bur[:, 1]
    nll_dataIO_pdfNO_Bur = arrIO_Bur[:, 3]
    nll_dataIO_pdfIO_Bur = arrIO_Bur[:, 5]

    arrNO_Gar = np.loadtxt("../results/Garching_10kpc_NO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modNO_Gar = arrNO_Gar[:, 0]
    tmaxNO_Gar = arrNO_Gar[:, 1]
    nll_dataNO_pdfNO_Gar = arrNO_Gar[:, 3]
    nll_dataNO_pdfIO_Gar = arrNO_Gar[:, 5]

    arrIO_Gar = np.loadtxt("../results/Garching_10kpc_IO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modIO_Gar = arrIO_Gar[:, 0]
    tmaxIO_Gar = arrIO_Gar[:, 1]
    nll_dataIO_pdfNO_Gar = arrIO_Gar[:, 3]
    nll_dataIO_pdfIO_Gar = arrIO_Gar[:, 5]

    colors0 = mpl.cm.gnuplot(np.linspace(0.2, 1, 8))
    colors1 = mpl.cm.viridis(np.linspace(0.2, 1, 8))

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    i = 0
    for x, y, m in zip(tmaxNO_Bur, nll_dataNO_pdfIO_Bur - nll_dataNO_pdfNO_Bur, modNO_Bur):
        if x == 0.01:
            ax[0].plot(x*1000, 2*y, "^", fillstyle="none", ms=8, color=colors0[i], label=f"{m}"+r"$M_\odot$")
        else:
            ax[0].plot(x*1000, 2*y, "^", fillstyle="none", ms=8, color=colors0[i])
        i += 1
        if i == 8:
            i = 0

    NmodsBur = int(len(tmaxIO_Bur) / 5)
    tarr = np.arange(10, 35, 5)
    #for g in range(5):
    #    ax[0].vlines(tarr[g]+1.5, np.min(2*(nll_dataNO_pdfIO_Bur - nll_dataNO_pdfNO_Bur)[g*NmodsBur:(g+1)*NmodsBur]), np.max(2*(nll_dataNO_pdfIO_Bur - nll_dataNO_pdfNO_Bur)[g*NmodsBur:(g+1)*NmodsBur]), linestyle="-", lw=2, color="blue")
    NmodsBur = int(len(tmaxIO_Gar) / 5)
    tarr = np.arange(10, 35, 5)
    #for g in range(5):
    #    ax[0].vlines(tarr[g]+2.0, np.min(2*(nll_dataNO_pdfIO_Gar - nll_dataNO_pdfNO_Gar)[g*NmodsBur:(g+1)*NmodsBur]), np.max(2*(nll_dataNO_pdfIO_Gar - nll_dataNO_pdfNO_Gar)[g*NmodsBur:(g+1)*NmodsBur]), linestyle="-", lw=2, color="crimson")

    i = 0
    for x, y, m in zip(tmaxNO_Gar, nll_dataNO_pdfIO_Gar - nll_dataNO_pdfNO_Gar, modNO_Gar):
        ax[0].plot(x*1000, 2*y, "v", fillstyle="none", ms=8, color=colors1[i])
        i += 1
        if i == 8:
            i = 0
    ax[0].set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=16)
    ax[0].set_ylabel(r"$\Delta \chi^2$", fontsize=16)
    ax[0].tick_params(axis="y", labelsize=15)
    ax[0].tick_params(axis="x", labelsize=15)
    ax[0].set_title("NO data", fontsize=16)
    ax[0].legend(prop={"size":10}, ncol=2, loc="lower right", frameon=True)
    ax[0].grid(True, linestyle=":")
    ax[0].set_ylim(0, 16)
    i = 0
    for x, y in zip(tmaxIO_Bur, nll_dataIO_pdfIO_Bur - nll_dataIO_pdfNO_Bur):
        ax[1].plot(x*1000, -2*y, "^", fillstyle="none", ms=8, color=colors0[i])
        i += 1
        if i == 8:
            i = 0
    i = 0
    for x, y, n in zip(tmaxIO_Gar, nll_dataIO_pdfIO_Gar - nll_dataIO_pdfNO_Gar, modNO_Gar):
        if x == 0.01:
            if n < 90000:
                ax[1].plot(x*1000, -2*y, "v", fillstyle="none", ms=8, color=colors1[i], label=f"LS220 EoS, {int((n-80000)/100.):.1f}"+r"$M_\odot$")
            else:
                ax[1].plot(x*1000, -2*y, "v", fillstyle="none", ms=8, color=colors1[i], label=f"Shen EoS, {int((n-80000)/100.):.1f}"+r"$M_\odot$")
        else:
            ax[1].plot(x*1000, -2*y, "v", fillstyle="none", ms=8, color=colors1[i])
        i += 1
        if i == 8:
            i = 0
    NmodsBur = int(len(tmaxIO_Bur) / 5)
    tarr = np.arange(10, 35, 5)
    #for g in range(5):
    #    ax[1].vlines(tarr[g]+1.5, np.min(2*(nll_dataIO_pdfNO_Bur - nll_dataIO_pdfIO_Bur)[g*NmodsBur:(g+1)*NmodsBur]), np.max(2*(nll_dataIO_pdfNO_Bur - nll_dataIO_pdfIO_Bur)[g*NmodsBur:(g+1)*NmodsBur]), linestyle="-", lw=2, color="blue")
    NmodsBur = int(len(tmaxIO_Gar) / 5)
    tarr = np.arange(10, 35, 5)
    #for g in range(5):
    #    ax[1].vlines(tarr[g]+2.0, np.min(2*(nll_dataIO_pdfNO_Gar - nll_dataIO_pdfIO_Gar)[g*NmodsBur:(g+1)*NmodsBur]), np.max(2*(nll_dataIO_pdfNO_Gar - nll_dataIO_pdfIO_Gar)[g*NmodsBur:(g+1)*NmodsBur]), linestyle="-", lw=2, color="crimson")

    ax[1].set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=16)
    ax[1].set_ylabel(r"$\Delta \chi^2$", fontsize=16)
    ax[1].tick_params(axis="y", labelsize=15)
    ax[1].tick_params(axis="x", labelsize=15)
    ax[1].set_title("IO data", fontsize=16)
    ax[1].legend(prop={"size":10}, ncol=2, loc="lower right", frameon=True)
    ax[1].grid(True, linestyle=":")
    ax[1].set_ylim(0, 16)

    plt.tight_layout()
    plt.savefig("../plots/compareModels_dark.pdf")
    plt.show()


if __name__ == "__main__":
    #drawTargetDistEthr3D()
    #drawDeltaTEthr()
    drawSensEthr()
    #drawSensTimeEthr()
    #drawSensTargetEthr()
    #drawSensDistEthr()
    #compareExps()
    #compareModels()

