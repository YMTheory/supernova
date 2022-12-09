import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def draw2D():
    arr = np.loadtxt("../data/CEvNS_quenched.txt")
    fig, ax = plt.subplots(figsize=(10, 4))
    #im = ax.imshow(arr, extent=[-20, 40, 0, 0.2], aspect="auto")
    #im = ax.imshow(arr, extent=[-20, 40, 0, 0.2], aspect="auto", norm=colors.LogNorm(vmin=arr.min(), vmax=arr.max()) )
    im = ax.imshow(arr, extent=[-20, 40, 0, 1.0], aspect="auto", cmap='viridis',  norm=colors.LogNorm(vmin=1e-6, vmax=arr.max()) )
    cb = plt.colorbar(im, ax = ax)
    cb.set_label(r"[MeV$^{-1}\cdot$s$^{-1}$]", fontsize=13)
    
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel(r"$E_\mathrm{vis}$ [MeV]", fontsize=14)
    ax.tick_params(axis="both", labelsize=13)
    
    ax.hlines(0.10, -20, 40, lw=2, color="red", linestyle="--")
    ax.hlines(0.15, -20, 40, lw=2, color="red", linestyle="-.")
    ax.set_ylim(0, 0.6)
    
    plt.tight_layout()
    plt.savefig("../plots/CEvNS_EvisT2D_quenched.pdf")
    plt.show()


def draw(Ethr):
    stepEvis = 0.01

    arr = np.loadtxt("../data/CEvNS_quenched.txt")

    startId = int(Ethr / stepEvis)
    endId = arr.shape[0] - startId
    print(Ethr, endId)
    
    arr1 = arr[0:endId, :]

    val = np.sum(arr1, axis=0) * stepEvis

    return val


def drawNuP(Emin, Emax):
    fig, ax = plt.subplots(figsize=(10, 4))
    tarr = np.arange(-30, 30, 1)
    for i in range(len(Emin)):
        E1, E2 = Emin[i], Emax[i]
        filename = f"./data/NuP_testInt_Emin{E1}MeVEmax{E2}MeV.txt"
        tmp_arr = np.loadtxt(filename)
        ax.plot(tarr, tmp_arr, lw=2, label=f"{E1}-{E2} MeV")
    ax.legend(prop={"size":12})
    ax.semilogy()
    plt.show()


def draw1D():
    c0 = draw(0.0)
    c1 = draw(0.10)
    c2 = draw(0.15)

    fig, ax = plt.subplots(figsize=(10, 4))
    t = np.arange(-20, 40, 1)
    ax.plot(t, c0, ":",  lw=2, label=r"E$_\mathrm{thr}$ = 0.00 MeV")
    ax.plot(t, c1, "-",  lw=2, label=r"E$_\mathrm{thr}$ = 0.10 MeV")
    ax.plot(t, c2, "--", lw=2, label=r"E$_\mathrm{thr}$ = 0.15 MeV")
    ax.legend(prop={"size":12})
    ax.grid(True, linestyle=":")
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("CEvNS-$^{12}$C counts per ms", fontsize=14)
    ax.tick_params(axis="both", labelsize=13)
    ax.semilogy()
    plt.tight_layout()
    plt.savefig("../plots/CEvNS_T1DEthr_quenched.pdf")
    plt.show()



def NuP_debug():
    arr = np.loadtxt("./data/NuP_check.txt")
    arr_low = arr[0, :]
    arr_hig = arr[1, :]
    arr_tot = arr[2, :]
    t = np.arange(-30, 30, 1)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, arr_low, lw=2, label="0.1-3 MeV")
    #ax.plot(t, arr_hig, lw=2, label="3-4 MeV")
    ax.plot(t, arr_tot, lw=2, label="0.1-4 MeV")
    ax.plot(t, arr_low+arr_hig, ":", lw=2, label="sum")
    ax.legend(prop={"size":12})
    #ax.semilogy()
    plt.show()







if __name__ == "__main__" :
    draw2D()
    #draw1D()
    #draw_combined()
    #drawNuP([0.1, 1, 3, 0.1], [1, 3, 5, 5])
    #NuP_debug()
