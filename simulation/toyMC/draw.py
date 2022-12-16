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

    asimovNO = 7.20
    asimovIO = 7.59
    
    _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", 0.15, 20, 10, end="tot_new_doFit_unbinnedNLL_PosiToyData")
    _, _, sens_dataNO_bin, _, _, sens_dataIO_bin    = read("Garching82703", "pESeESIBD", 0.15, 20, 10, end="tot_new_doFit_binnedNLL_PosiToyData")

    fig, ax = plt.subplots(figsize=(9, 4))
    h1 = ax.hist(sens_dataNO,   bins=100, range=(-20, 20),   alpha=0.5, label="true NO (unbinned)", color="blue",)
    h2 = ax.hist(-sens_dataIO,  bins=100, range=(-20, 20),  alpha=0.5, label="true IO (unbinned)", color="red")
    h3 = ax.hist(sens_dataNO_bin,   bins=100, range=(-20, 20),   alpha=0.5, label="true NO (0.01 ms binning)", histtype="step", lw=2, color="green",)
    h4 = ax.hist(-sens_dataIO_bin,  bins=100, range=(-20, 20),  alpha=0.5, label="true IO (0.01 ms binning)", histtype="step", lw=2, color="orange")
    ax.set_xlabel(r"$\Delta\chi^2$", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=0)
    #ax.set_ylim(0, 7000)

    l1 = ax.vlines(-np.median(sens_dataIO), 0, 2000, linestyle=":", lw=2.5, color="red")
    l2 = ax.vlines(-np.median(sens_dataIO_bin), 0, 2500, linestyle="-.", lw=2.5, color="orange")
    l3 = ax.vlines(-asimovIO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(-np.median(sens_dataIO)-4, 1000, f"{-np.median(sens_dataIO):.1f}", fontsize=15, color="red")
    ax.text(-np.median(sens_dataIO_bin)-4, 1300, f"{-np.median(sens_dataIO_bin):.1f}", fontsize=15, color="orange")
    ax.text(-np.median(sens_dataIO)+2, 1500, f"{-asimovIO}", fontsize=15, color="black")

    ax.vlines(np.median(sens_dataNO),  0, 2000, linestyle=":", lw=2.5, color="blue")
    ax.vlines(np.median(sens_dataNO_bin),  0, 50, linestyle="-.", lw=2.5, color="green")
    ax.vlines(asimovNO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(np.median(sens_dataNO)-2, 1000, f"{np.median(sens_dataNO):.1f}", fontsize=15, color="blue")
    ax.text(np.median(sens_dataNO_bin)-4, 1300, f"{np.median(sens_dataNO_bin):.1f}", fontsize=15, color="green")
    ax.text(np.median(sens_dataNO)+2, 1500, f"{asimovNO}", fontsize=15, color="black")
    
    #lg1 = ax.legend([l1, l3], ["median value of Monte Carlo", "Asimov dataset"], loc="upper right", frameon=True)
    lg1 = ax.legend([l1, l2, l3], ["unbinned", "0.01 ms binned", "Asimov"], loc="upper right", frameon=True)
    plt.legend(prop={"size":12}, loc="upper left", frameon=True, ncol=2)
    plt.gca().add_artist(lg1)

    plt.tight_layout()
    #plt.savefig("Garching82703_pESeESIBD_10kpc_noFit_scale10.pdf")
    plt.show()





if __name__ == "__main__":
    draw()





