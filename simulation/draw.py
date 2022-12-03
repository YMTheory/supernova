import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def draw2D():
    arr = np.loadtxt("CEvNS_quenched.txt")
    fig, ax = plt.subplots(figsize=(8, 4))
    #im = ax.imshow(arr, extent=[-20, 40, 0, 0.2], aspect="auto")
    #im = ax.imshow(arr, extent=[-20, 40, 0, 0.2], aspect="auto", norm=colors.LogNorm(vmin=arr.min(), vmax=arr.max()) )
    im = ax.imshow(arr, extent=[-20, 40, 0, 0.5], aspect="auto", norm=colors.LogNorm(vmin=1e-12, vmax=arr.max()) )
    cb = plt.colorbar(im, ax = ax)
    cb.set_label(r"[MeV$^{-1}\cdot$s$^{-1}$]", fontsize=13)
    
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel(r"$E_\mathrm{vis}$ [MeV]", fontsize=14)
    ax.tick_params(axis="both", labelsize=13)
    
    ax.hlines(0.10, -20, 40, lw=2, color="red", linestyle="--")
    ax.hlines(0.15, -20, 40, lw=2, color="red", linestyle="-.")
    
    plt.tight_layout()
    plt.savefig("CEvNS_EvisT2D_quenched.pdf")
    plt.show()


def draw(Ethr):
    stepEvis = 0.005

    arr = np.loadtxt("NuP_quenched.txt")

    startId = int(Ethr / stepEvis)
    endId = arr.shape[0] - startId
    print(Ethr, endId)
    
    arr1 = arr[0:endId, :]

    val = np.sum(arr1, axis=0) * stepEvis

    return val




def draw1D():
    c0 = draw(0.0)
    c1 = draw(0.10)
    c2 = draw(0.15)

    fig, ax = plt.subplots(figsize=(8, 4))
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
    plt.savefig("NuP_T1DEthr_quenched.pdf")
    plt.show()


def draw_combined():
    stepEvis = 0.005

    arr = np.loadtxt("CEvNS_quenched.txt")
    c0 = draw(0.0)
    c1 = draw(0.10)
    c2 = draw(0.15)

    fig, ax = plt.subplots(2, 1, figsize=(8, 6))

    im = ax[1].imshow(arr, extent=[-20, 40, 0, 0.5], aspect="auto", norm=colors.LogNorm(vmin=1e-9, vmax=arr.max()) )
    cb = plt.colorbar(im, ax = ax[1])
    cb.set_label(r"[MeV$^{-1}\cdot$s$^{-1}$]", fontsize=13)
    ax[1].set_xlabel("post-bounce time [ms]", fontsize=14)
    ax[1].set_ylabel(r"$E_\mathrm{vis}$ [MeV]", fontsize=14)
    ax[1].tick_params(axis="both", labelsize=13)
    ax[1].hlines(0.10, -20, 40, lw=2, color="red", linestyle="--")
    ax[1].hlines(0.15, -20, 40, lw=2, color="red", linestyle="-.")

    t = np.arange(-20, 40, 1)
    ax[0].plot(t, c0, ":",  lw=2, label=r"E$_\mathrm{thr}$ = 0.00 MeV")
    ax[0].plot(t, c1, "-",  lw=2, label=r"E$_\mathrm{thr}$ = 0.10 MeV")
    ax[0].plot(t, c2, "--", lw=2, label=r"E$_\mathrm{thr}$ = 0.15 MeV")
    ax[0].legend(prop={"size":12})
    ax[0].grid(True, linestyle=":")
    #ax[0, 1].set_xlabel("post-bounce time [ms]", fontsize=14)
    ax[0].set_ylabel("CEvNS-$^{12}$C counts per ms", fontsize=14)
    ax[0].tick_params(axis="both", labelsize=13)
    ax[0].semilogy()

    #plt.subplots_adjust(hspace=0.01)
    plt.tight_layout()
    plt.savefig("CEvNS_quenched.pdf")
    plt.show()



if __name__ == "__main__" :
    #draw2D()
    draw1D()
    #draw_combined()

