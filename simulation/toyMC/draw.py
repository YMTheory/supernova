import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import ROOT
import uproot as up


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
    asimovNO = 7.3
    asimovIO = 8.4

    Et = 0.10
    fig, ax = plt.subplots(figsize=(9, 4))
    _, _, sens_dataNO, _, _, sens_dataIO                = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    h1 = ax.hist(sens_dataNO,   bins=100, range=(-30, 30),   alpha=0.5, label="true NO", color="blue",)
    h2 = ax.hist(-sens_dataIO,  bins=100, range=(-30, 30),  alpha=0.5, label="true IO", color="red")


    conts1, edges1 = np.histogram(sens_dataNO, bins=100, range=(-30, 30))
    conts2, edges2 = np.histogram(-sens_dataIO, bins=100, range=(-30, 30))
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
    ax.text(16, 2000, "gauss", fontsize=15)



    
    ax.set_xlabel(r"$\Delta\chi^2$", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=0)
    #ax.set_ylim(0, 7000)

    l1 = ax.vlines(-np.median(sens_dataIO), 0, 3500, linestyle=":", lw=2.5, color="red")
    l3 = ax.vlines(-asimovIO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(-np.median(sens_dataIO)-4, 1000, f"{-np.median(sens_dataIO):.1f}", fontsize=15, color="red")
    ax.text(-asimovIO+2, 1500, f"{-asimovIO:.1f}", fontsize=15, color="black")

    ax.vlines(np.median(sens_dataNO),  0, 3500, linestyle=":", lw=2.5, color="blue")
    ax.vlines(asimovNO, 0, 2000, linestyle="--", lw=2.5, color="black")
    ax.text(np.median(sens_dataNO)+2, 1000, f"{np.median(sens_dataNO):.1f}", fontsize=15, color="blue")
    ax.text(asimovNO+2, 1500, f"{asimovNO:.1f}", fontsize=15, color="black")
    ax.set_ylim(0, 6000)
    
    lg1 = ax.legend([l1, l3], ["median value of Monte Carlo", "Asimov dataset"], loc="upper right", frameon=True)
    plt.legend(prop={"size":12}, loc="upper left", frameon=True, ncol=3)
    plt.gca().add_artist(lg1)

    plt.tight_layout()
    plt.savefig(f"./plots/Garching82703_pESeESIBD_10kpc_0.1MeV_noC14bkg.png")
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

    fig, ax = plt.subplots(figsize=(9, 4))
    #fig, (ax1, ax) = plt.subplots(2, 1, figsize=(9, 5))

    Et_arr = [0.00, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, ]
    sens_NO_arr, sens_IO_arr = [], []
    for Et in Et_arr:
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", Et, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
        sens_NO_arr.append(np.median(sens_dataNO))
        sens_IO_arr.append(np.median(sens_dataIO))

    ax.plot(Et_arr, sens_NO_arr, "^-", ms=8, lw=2, label="NO: w/o C14")
    ax.plot(Et_arr, sens_IO_arr, "v-", ms=8, lw=2, label="IO: w/o C14")

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
    ax.fill_between(Et_arr, sens_IO_arrl, sens_IO_arrh, color="orange", alpha=0.5, label="IO: w/ C14")


    #_, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBDCEvNS", 0.00, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #ax.plot(0, np.median(sens_dataNO), "v", ms=8, color="darkviolet")
    #ax.plot(0, np.median(sens_dataIO), "^", ms=8, color="darkviolet")
    #_, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBD", 0.10, 20, 10, end="tot_nodoFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #ax.hlines(np.median(sens_dataNO), 0, 0.2, linestyle="--", lw=2, alpha=0.6, color="green")
    #ax.hlines(np.median(sens_dataIO), 0, 0.2, linestyle="-.", lw=2, alpha=0.6, color="green")


    #ax.semilogx()
    ax.set_xlabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax.set_ylabel(r"$\Delta\chi^2$", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
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
    plt.savefig(f"./plots/Garching82703_pESeESIBD+CEvNS_10kpc_C14bkg.png")
    plt.show()

def draw_models():

    fig, ax = plt.subplots(figsize=(9, 5))
    colors0 = mpl.cm.gnuplot(np.linspace(0.1,1, 8))
    sens_NO_arr, sens_IO_arr = [], []
    mods = ["81123", "82503", "82703", "84003", "91123", "92503", "92703", "94003"]
    labels = [r"Shen EoS, $11.2M_\odot$", r"Shen EoS, $25M_\odot$", r"Shen EoS, $27M_\odot$", r"Shen EoS, $40M_\odot$", r"LS200 EoS, $11.2M_\odot$", r"LS200 EoS, $25M_\odot$", r"LS200 EoS, $27M_\odot$", r"LS200 EoS, $40M_\odot$", ]
    tmax_arr = [15, 20, 25, 30]
    for i, mod in enumerate(mods):
        m_sens_dataNO, m_sens_dataIO = [], []
        for tmax in tmax_arr:
            _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching"+mod, "pESeESIBD", 0.10, tmax, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
            m_sens_dataNO.append(np.median(sens_dataNO))
            m_sens_dataIO.append(np.median(sens_dataIO))
            
        ax.plot(tmax_arr, m_sens_dataNO, "^", color=colors0[i], ms=8, lw=2, label=labels[i])
        ax.plot(tmax_arr, m_sens_dataIO, "v", color=colors0[i], ms=8, lw=2, fillstyle="none")

    #_, _, sens_dataNO, _, _, sens_dataIO            = read("Garching82703", "pESeESIBDCEvNS", 0.00, 20, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
    #ax.plot(0, np.median(sens_dataNO), "v", ms=8, color="red")
    #ax.plot(0, np.median(sens_dataIO), "^", ms=8, color="red")



    sens_NO_arr, sens_IO_arr = [], []
    sens_NO_arr.append(0)
    sens_IO_arr.append(0)
    mods = ["82503", "82703", "84003", "91123", "92503", "92703", "94003"]
    for i, mod in enumerate(mods):
        m_sens_dataNO, m_sens_dataIO = [], []
        _, _, sens_dataNO, _, _, sens_dataIO            = read("Garching"+mod, "pESeES", 0.10, 10, 10, end="tot_doFit_binnedData_unbinnedNLL_noC14bkg_PosiToyData")
        m_sens_dataNO.append(np.median(sens_dataNO))
        m_sens_dataIO.append(np.median(sens_dataIO))
            
        ax.plot([10], m_sens_dataNO, "o", color=colors0[i], ms=8, lw=2)
        ax.plot([10], m_sens_dataIO, "s", color=colors0[i], ms=8, lw=2, fillstyle="none")




    #ax.semilogx()
    ax.set_xlabel(r"$t_\mathrm{max}$ [ms]", fontsize=16)
    ax.set_ylabel(r"$\Delta\chi^2$", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.grid(True, linestyle=":")
    ax.legend(prop={"size":12}, loc="upper left", frameon=True, ncol=2)

    ax.text(9.5, 4.5, "Normal ordering", fontsize=15)
    ax.text(9.5, 6.0, "Inveted ordering", fontsize=15)

    plt.tight_layout()
    plt.savefig(f"./plots/Garchings_pESeESIBD_202530ms_10kpc_noC14bkg.png")
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
    cb.set_label(r"$\Delta \chi^2$", fontsize=14)
    im1 = ax1.imshow(m_sens_dataIO, extent=[4.5, 10.5, 0.095, 0.155], aspect="auto")
    cb1 = plt.colorbar(im1, ax=ax1)
    cs = ax.contour(X, Y, m_sens_dataNO, levels=[9, 16,  25], colors=[ "red", "darkorange", "darkviolet"])
    cs1 = ax1.contour(X, Y, m_sens_dataIO ,levels=[9, 16, 25],colors=[ "red", "darkorange", "darkviolet"] )


    cb1.set_label(r"$\Delta \chi^2$", fontsize=14)

    #ax.semilogx()
    ax.set_xlabel("Distance [kpc]", fontsize=16)
    ax.set_ylabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax1.set_xlabel("Distance [kpc]", fontsize=16)
    ax1.set_ylabel(r"$E_\mathrm{thr}$ [MeV]", fontsize=16)
    ax1.tick_params(axis="x", labelsize=15)
    ax1.tick_params(axis="y", labelsize=15)


    plt.tight_layout()
    plt.savefig(f"./plots/Garching82703_pESeESIBD_20ms_distances_C14bkglow.png")
    plt.show()


if __name__ == "__main__":
    #draw()
    drawC14()
    #draw_Ethr()
    #drawNevent()
    #draw_Tmax()
    #draw_models()
    #draw_distance()




