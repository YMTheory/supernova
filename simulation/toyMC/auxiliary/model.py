import numpy as np
import matplotlib.pyplot as plt
import uproot as up
from matplotlib import pyplot as plt, cm
from matplotlib import colors
import matplotlib as mpl
plt.style.use("science")



if __name__ == "__main__":

    arrNO_Bur = np.loadtxt("../results/Burrows_10kpc_NO_noC14_pESeES_scale5_tmax15202530ms_MOSensAsimov_JUNO.csv")
    #arrNO_Bur = np.loadtxt("../results/Burrows_10kpc_NO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modNO_Bur = arrNO_Bur[:, 0]
    tmaxNO_Bur = arrNO_Bur[:, 1]
    nll_dataNO_pdfNO_Bur = arrNO_Bur[:, 3]
    nll_dataNO_pdfIO_Bur = arrNO_Bur[:, 5]

    arrIO_Bur = np.loadtxt("../results/Burrows_10kpc_IO_noC14_pESeES_scale5_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modIO_Bur = arrIO_Bur[:, 0]
    tmaxIO_Bur = arrIO_Bur[:, 1]
    nll_dataIO_pdfNO_Bur = arrIO_Bur[:, 3]
    nll_dataIO_pdfIO_Bur = arrIO_Bur[:, 5]

    arrNO_Gar = np.loadtxt("../results/Garching_10kpc_NO_noC14_pESeES_scale5_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modNO_Gar = arrNO_Gar[:, 0]
    tmaxNO_Gar = arrNO_Gar[:, 1]
    nll_dataNO_pdfNO_Gar = arrNO_Gar[:, 3]
    nll_dataNO_pdfIO_Gar = arrNO_Gar[:, 5]

    arrIO_Gar = np.loadtxt("../results/Garching_10kpc_IO_noC14_pESeES_scale5_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modIO_Gar = arrIO_Gar[:, 0]
    tmaxIO_Gar = arrIO_Gar[:, 1]
    nll_dataIO_pdfNO_Gar = arrIO_Gar[:, 3]
    nll_dataIO_pdfIO_Gar = arrIO_Gar[:, 5]

    arrNO_Bur0 = np.loadtxt("../results/Burrows_10kpc_NO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modNO_Bur0 = arrNO_Bur0[:, 0]
    tmaxNO_Bur0 = arrNO_Bur0[:, 1]
    nll_dataNO_pdfNO_Bur0 = arrNO_Bur0[:, 3]
    nll_dataNO_pdfIO_Bur0 = arrNO_Bur0[:, 5]

    arrIO_Bur0 = np.loadtxt("../results/Burrows_10kpc_IO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modIO_Bur0 = arrIO_Bur0[:, 0]
    tmaxIO_Bur0 = arrIO_Bur0[:, 1]
    nll_dataIO_pdfNO_Bur0 = arrIO_Bur0[:, 3]
    nll_dataIO_pdfIO_Bur0 = arrIO_Bur0[:, 5]

    arrNO_Gar0 = np.loadtxt("../results/Garching_10kpc_NO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modNO_Gar0 = arrNO_Gar0[:, 0]
    tmaxNO_Gar0 = arrNO_Gar0[:, 1]
    nll_dataNO_pdfNO_Gar0 = arrNO_Gar0[:, 3]
    nll_dataNO_pdfIO_Gar0 = arrNO_Gar0[:, 5]

    arrIO_Gar0 = np.loadtxt("../results/Garching_10kpc_IO_noC14_pESeESIBD_tmax15202530ms_MOSensAsimov_JUNO.csv")
    modIO_Gar0 = arrIO_Gar0[:, 0]
    tmaxIO_Gar0 = arrIO_Gar0[:, 1]
    nll_dataIO_pdfNO_Gar0 = arrIO_Gar0[:, 3]
    nll_dataIO_pdfIO_Gar0 = arrIO_Gar0[:, 5]



    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    pos = [10, 9, 8, 7, 6, 5]
    dchi2_dataNO_Bur = 2 * (nll_dataNO_pdfIO_Bur - nll_dataNO_pdfNO_Bur)
    dchi2_dataNO_Bur = dchi2_dataNO_Bur.reshape((6, 8))
    dchi2_dataNO_Bur = dchi2_dataNO_Bur.T
    violin_parts = ax[0].violinplot(dataset=dchi2_dataNO_Bur , showmeans=True, positions=pos, widths=1)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('royalblue')
        pc.set_edgecolor('black')
    dchi2_dataNO_Gar = 2 * (nll_dataNO_pdfIO_Gar - nll_dataNO_pdfNO_Gar)
    dchi2_dataNO_Gar = dchi2_dataNO_Gar.reshape((6, 8))
    dchi2_dataNO_Gar = dchi2_dataNO_Gar.T
    violin_parts = ax[0].violinplot(dataset=dchi2_dataNO_Gar, showmeans=True, positions=pos, widths=1)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('green')
        pc.set_edgecolor('black')
    pos0 = [10, 15, 20, 25, 30]
    dchi2_dataNO_Bur0 = 2 * (nll_dataNO_pdfIO_Bur0 - nll_dataNO_pdfNO_Bur0)
    dchi2_dataNO_Bur0 = dchi2_dataNO_Bur0.reshape((5, 8))
    dchi2_dataNO_Bur0 = dchi2_dataNO_Bur0.T
    ax[0].violinplot(dataset=dchi2_dataNO_Bur0, showmeans=True, positions=pos0, widths=2)
    dchi2_dataNO_Gar0 = 2 * (nll_dataNO_pdfIO_Gar0 - nll_dataNO_pdfNO_Gar0)
    dchi2_dataNO_Gar0 = dchi2_dataNO_Gar0.reshape((5, 8))
    dchi2_dataNO_Gar0 = dchi2_dataNO_Gar0.T
    ax[0].violinplot(dataset=dchi2_dataNO_Gar0, showmeans=True, positions=pos0, widths=2)
    #ax[0].set_xticks([0, 1, 2, 3, 4], [10, 15, 20, 25, 30])
    ax[0].set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=16)
    ax[0].set_ylabel(r"$\Delta\chi^2$", fontsize=16)
    ax[0].tick_params(axis="both", labelsize=15)
    ax[0].set_ylim(0, 15)
    ax[0].grid(True, linestyle=":", alpha=0.5)
    ax[0].set_title("NO data", fontsize=16)

    dchi2_dataIO_Bur = -2 * (nll_dataIO_pdfIO_Bur - nll_dataIO_pdfNO_Bur)
    dchi2_dataIO_Bur = dchi2_dataIO_Bur.reshape((6, 8))
    dchi2_dataIO_Bur = dchi2_dataIO_Bur.T
    violin_parts = ax[1].violinplot(dataset=dchi2_dataIO_Bur , showmeans=True, positions=pos, widths=1)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('royalblue')
        pc.set_edgecolor('black')
    dchi2_dataIO_Gar = -2 * (nll_dataIO_pdfIO_Gar - nll_dataIO_pdfNO_Gar)
    dchi2_dataIO_Gar = dchi2_dataIO_Gar.reshape((6, 8))
    dchi2_dataIO_Gar = dchi2_dataIO_Gar.T
    violin_parts = ax[1].violinplot(dataset=dchi2_dataIO_Gar , showmeans=True, positions=pos, widths=1)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('green')
        pc.set_edgecolor('black')
    dchi2_dataIO_Bur0 = -2 * (nll_dataIO_pdfIO_Bur0 - nll_dataIO_pdfNO_Bur0)
    dchi2_dataIO_Bur0 = dchi2_dataIO_Bur0.reshape((5, 8))
    dchi2_dataIO_Bur0 = dchi2_dataIO_Bur0.T
    ax[1].violinplot(dataset=dchi2_dataIO_Bur0 , showmeans=True, positions=pos0, widths=2)
    dchi2_dataIO_Gar0 = -2 * (nll_dataIO_pdfIO_Gar0 - nll_dataIO_pdfNO_Gar0)
    dchi2_dataIO_Gar0 = dchi2_dataIO_Gar0.reshape((5, 8))
    dchi2_dataIO_Gar0 = dchi2_dataIO_Gar0.T
    ax[1].violinplot(dataset=dchi2_dataIO_Gar0 , showmeans=True, positions=pos0, widths=2)
    ax[1].set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=16)
    ax[1].set_ylabel(r"$\Delta\chi^2$", fontsize=16)
    ax[1].tick_params(axis="both", labelsize=15)
    ax[1].set_ylim(0, 15)
    ax[1].grid(True, linestyle=":", alpha=0.5)
    ax[1].plot(10, 100, "-", lw=3, color="royalblue", label="Princeton (ES only)")
    ax[1].plot(10, 100, "-", lw=3, color="green", label="Garching (ES only)")
    ax[1].plot(10, 100, "-", lw=3, color="orange", label="Princeton (default)")
    ax[1].plot(10, 100, "-", lw=3, color="red", label="Garching (default)")
    ax[1].legend(prop={"size":12}, frameon=True, loc="lower right", ncol=2)
    ax[1].set_title("IO data", fontsize=16)

    plt.tight_layout()
    plt.savefig("../plots/modelDepend_ESonlyScale5.pdf")
    plt.show()





