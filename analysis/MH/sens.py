import numpy as np
import matplotlib.pyplot as plt
plt.style.use('science')


def read_data(filename):
    model, dchi2, N1, N2, N3 = [], [], [], [], []
    with open(filename, "r") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Ndata = len(data)
            model.append(int(data[0]))
            dchi2.append(float(data[2]))
            if Ndata >= 4:
                N1.append(float(data[3]))
            if Ndata >= 5:
                N2.append(float(data[4]))
            if Ndata >= 6:
                N3.append(float(data[5]))
    dchi2 = np.array(dchi2)
    return model, dchi2, N1, N2, N3


import sys
if __name__ == "__main__" :

    Ethr = 0.20
    mod = None
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv), 1):
            if sys.argv[i] == "-Ethr":
                Ethr = float(sys.argv[i+1])
            if sys.argv[i] == "-model" :
                mod = sys.argv[i+1]

    filename = "/junofs/users/miaoyu/supernova/analysis/MH/results/res_%s_NO_pESIBDeES_%.2fMeV.txt"%(mod, Ethr)
    modelNO, dchi2NO, N1NO, N2NO, N3NO = read_data(filename)
    filename = "/junofs/users/miaoyu/supernova/analysis/MH/results/res_%s_IO_pESIBDeES_%.2fMeV.txt" %(mod, Ethr)
    modelIO, dchi2IO, N1IO, N2IO, N3IO = read_data(filename)

    print("Maximum, medium and minimum sensitivity for NO data: %.2f, %.2f, %.2f"%(np.max(dchi2NO), np.median(dchi2NO), np.min(dchi2NO)))
    print("Maximum, medium and minimum sensitivity for IO data: %.2f, %.2f, %.2f"%(np.max(dchi2IO), np.median(dchi2IO), np.min(dchi2IO)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    stitle = "Garching-10kpc-pESeESIBD-%.2fMeV"%(Ethr)
    st = fig.suptitle(stitle, fontsize=16)
    if mod == "Garching32":
        E = range(0, 32, 1)
    elif mod == "Burrows2D6000":
        E = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26]
    ax1.plot(E, dchi2NO, "o-", ms=6, lw=2, label="NO")
    ax1.plot(E, dchi2IO, "s-", ms=6, lw=2, label="IO")

    ax2.plot(E, N1NO, "X-", ms=6, lw=2, label="pES")
    ax2.plot(E, N2NO, "^-", color="seagreen", ms=6, lw=2, label="IBD, NO")
    ax2.plot(E, N3NO, "v-", color="orange", ms=6, lw=2, label="eES, NO")
    ax2.plot(E, N2IO, "^--", color="seagreen", fillstyle="none", ms=6, lw=2, label="IBD, IO")
    ax2.plot(E, N3IO, "v--", color="orange",   fillstyle="none", ms=6, lw=2, label="eES, IO")

    if mod == "Garching32":
        ax1.set_xlabel("Model No.", fontsize=15)
    elif mod == "Burrows2D6000":
        ax1.set_xlabel("Solar masses", fontsize=15)
    ax1.set_ylabel(r"$\Delta \chi^2$", fontsize=15)
    ax1.legend(prop={"size":16})
    ax1.tick_params(axis='both', which='major', labelsize=16, labelcolor="black")

    ax2.set_xlabel("Model No.", fontsize=15)
    ax2.set_ylabel("Event number in ROI", fontsize=15)
    ax2.legend(loc="upper left", ncol=3,  prop={"size":16})
    ax2.set_ylim(0, 40)
    ax2.tick_params(axis='both', which='major', labelsize=16, labelcolor="black")

    plt.tight_layout()
    plt.savefig("/junofs/users/miaoyu/supernova/analysis/MH/results/%s-10kpc-pESeESIBD-%.2fMeV.pdf"%(mod, Ethr))
    plt.show()




