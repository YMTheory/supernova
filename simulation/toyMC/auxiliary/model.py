import numpy as np
import matplotlib.pyplot as plt
import uproot as up
from matplotlib import pyplot as plt, cm
from matplotlib import colors
import matplotlib as mpl
plt.style.use("science")



if __name__ == "__main__":

    arrNO_Bur = np.loadtxt("../results/Burrows_10kpc_NO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modNO_Bur = arrNO_Bur[:, 0]
    tmaxNO_Bur = arrNO_Bur[:, 1]
    nll_dataNO_pdfNO_Bur = arrNO_Bur[:, 3]
    nll_dataNO_pdfIO_Bur = arrNO_Bur[:, 5]

    arrIO_Bur = np.loadtxt("../results/Burrows_10kpc_IO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
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

    colors0 = mpl.cm.gnuplot(np.linspace(0.1, 1, 8))
    colors1 = mpl.cm.viridis(np.linspace(0.1, 1, 8))
    #pl0 = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)
    #pl1 = sns.color_palette("dark:#5A9_r", as_cmap=True)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    dchi2_dataNO_Bur = 2 * (nll_dataNO_pdfIO_Bur - nll_dataNO_pdfNO_Bur)
    dchi2_dataNO_Bur = dchi2_dataNO_Bur.reshape((5, 8))
    dchi2_dataNO_Bur = dchi2_dataNO_Bur.T
    ax[0].violinplot(dataset=dchi2_dataNO_Bur , showmeans=True, positions=[10, 15, 20, 25, 30], widths=3)
    dchi2_dataNO_Gar = 2 * (nll_dataNO_pdfIO_Gar - nll_dataNO_pdfNO_Gar)
    dchi2_dataNO_Gar = dchi2_dataNO_Gar.reshape((5, 8))
    dchi2_dataNO_Gar = dchi2_dataNO_Gar.T
    ax[0].violinplot(dataset=dchi2_dataNO_Gar, showmeans=True, positions=[10, 15, 20, 25, 30], widths=3)
    #ax[0].set_xticks([0, 1, 2, 3, 4], [10, 15, 20, 25, 30])
    ax[0].set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=16)
    ax[0].set_ylabel(r"$\Delta\chi^2$", fontsize=16)
    ax[0].tick_params(axis="both", labelsize=15)
    ax[0].set_ylim(0, 15)
    ax[0].plot(0, 100, "-", lw=3, color="royalblue", label="Burrows et al.")
    ax[0].plot(0, 100, "-", lw=3, color="green", label="Garching group")
    ax[0].legend(prop={"size":14}, frameon=True, loc="lower right")
    ax[0].grid(True, linestyle=":", alpha=0.5)

    dchi2_dataIO_Bur = -2 * (nll_dataIO_pdfIO_Bur - nll_dataIO_pdfNO_Bur)
    dchi2_dataIO_Bur = dchi2_dataIO_Bur.reshape((5, 8))
    dchi2_dataIO_Bur = dchi2_dataIO_Bur.T
    ax[1].violinplot(dataset=dchi2_dataIO_Bur , showmeans=True, positions=[10, 15, 20, 25, 30], widths=3)
    dchi2_dataIO_Gar = -2 * (nll_dataIO_pdfIO_Gar - nll_dataIO_pdfNO_Gar)
    dchi2_dataIO_Gar = dchi2_dataIO_Gar.reshape((5, 8))
    dchi2_dataIO_Gar = dchi2_dataIO_Gar.T
    ax[1].violinplot(dataset=dchi2_dataIO_Gar , showmeans=True, positions=[10, 15, 20, 25, 30], widths=3)
    ax[1].set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=16)
    ax[1].set_ylabel(r"$\Delta\chi^2$", fontsize=16)
    ax[1].tick_params(axis="both", labelsize=15)
    ax[1].set_ylim(0, 15)
    ax[1].grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig("../plots/modelDepend_pESeESIBD.pdf")
    plt.show()





