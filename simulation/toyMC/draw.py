import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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




def draw():
    # Ethr = 0.15 MeV
    asimovNO = 7.20
    asimovIO = 7.59
    # Ethr = 0.10 MeV
    asimovNO = 7.7
    asimovIO = 8.7

    Et = 0.10
    
    _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_PosiToyData")
    _, _, sens_dataNO_low, _, _, sens_dataIO_low    = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_C14bkglow_PosiToyData")
    _, _, sens_dataNO_high, _, _, sens_dataIO_high    = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_C14bkghigh_PosiToyData")

    fig, ax = plt.subplots(figsize=(9, 4))
    h1 = ax.hist(sens_dataNO,   bins=100, range=(-30, 30),   alpha=0.5, label="true NO (no C14)", color="blue",)
    h2 = ax.hist(-sens_dataIO,  bins=100, range=(-30, 30),  alpha=0.5, label="true IO (no C14)", color="red")
    h3 = ax.hist(sens_dataNO_low,   bins=100, range=(-30, 30),   alpha=0.5, label="true NO (low C14)", histtype="step", lw=2, color="green",)
    h4 = ax.hist(-sens_dataIO_low,  bins=100, range=(-30, 30),  alpha=0.5, label="true IO (low C14)", histtype="step", lw=2, color="orange")
    h5 = ax.hist(sens_dataNO_high,   bins=100, range=(-30, 30),  linestyle=":",  alpha=0.5, label="true NO (high C14)", histtype="step", lw=2, color="darkviolet",)
    h6 = ax.hist(-sens_dataIO_high,  bins=100, range=(-30, 30),  linestyle=":", alpha=0.5, label="true IO (high C14)", histtype="step", lw=2, color="brown")
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


def draw_Ethr():

    fig, ax = plt.subplots(figsize=(9, 4))
    sens_NO_arr, sens_IO_arr = [], []
    Et_arr = [0.05, 0.08, 0.10, 0.12, 0.15]
    for Et in Et_arr:
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
        sens_NO_arr.append(np.median(sens_dataNO))
        sens_IO_arr.append(np.median(sens_dataIO))

    ax.plot(Et_arr, sens_NO_arr, "v-", ms=8, lw=2, label="Normal ordering")
    ax.plot(Et_arr, sens_IO_arr, "^-", ms=8, lw=2, label="Inverted ordering")

    ax.set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax.set_ylabel(r"$\Delta\chi^2$", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.grid(True, linestyle=":")
    ax.legend(prop={"size":12}, loc="upper right", frameon=True)

    plt.tight_layout()
    #plt.savefig(f"./plots/Garching82703_pESeESIBD_10kpc_{Et:.2f}MeV_withC14bkg.png")
    plt.show()


if __name__ == "__main__":
    #draw()
    draw_Ethr()





