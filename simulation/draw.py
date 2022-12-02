import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def draw2D():
    arr = np.loadtxt("CEvNS.txt")
    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(arr, extent=[-20, 40, 0, 0.2], aspect="auto", norm=colors.LogNorm(vmin=arr.min(), vmax=arr.max()) )
    cb = plt.colorbar(im, ax = ax)
    cb.set_label(r"[MeV$^{-1}\cdot$s$^{-1}$]", fontsize=13)
    
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel(r"$E_\mathrm{vis}$ [MeV]", fontsize=14)
    ax.tick_params(axis="both", labelsize=13)
    
    ax.hlines(0.10, -20, 40, lw=2, color="red", linestyle="--")
    ax.hlines(0.15, -20, 40, lw=2, color="red", linestyle="-.")
    
    plt.tight_layout()
    plt.savefig("CEvNS_EvisT2D.pdf")
    plt.show()


def draw(Ethr):
    stepEvis = 0.005

    arr = np.loadtxt("CEvNS.txt")

    startId = int(Ethr / stepEvis)
    endId = arr.shape[0] - startId
    
    arr1 = arr[0:endId, :]

    val = np.sum(arr1, axis=0) * stepEvis

    return val




if __name__ == "__main__":
    c1 = draw(0.10)
    c2 = draw(0.15)

    fig, ax = plt.subplots(figsize=(6, 4))
    t = np.arange(-20, 40, 1)
    ax.plot(t, c1, "-",  lw=2, label=r"E$_\mathrm{thr}$ = 0.10 MeV")
    ax.plot(t, c2, "--", lw=2, label=r"E$_\mathrm{thr}$ = 0.15 MeV")
    ax.legend(prop={"size":12})
    ax.grid(True, linestyle=":")
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("CEvNS-$^{12}$C counts per ms", fontsize=14)
    ax.tick_params(axis="both", labelsize=13)
    plt.tight_layout()
    plt.savefig("CEvNS_T1DEthr.pdf")
    plt.show()
