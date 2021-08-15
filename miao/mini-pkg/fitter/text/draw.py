import numpy as np
import matplotlib.pyplot as plt

def read_nll2T(filename):
    deltaT, nsig, nll = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            nll.append(float(data[0]))
            nsig.append(float(data[1]))
            deltaT.append(float(data[2]))

    deltaT = np.array(deltaT)
    nsig = np.array(nsig)
    nll = 2 * np.array(nll)

    return deltaT, nsig, nll





if __name__ == "__main__":
    deltaT1D, nsig1D, nll1D = read_nll2T("fineNll2T_dataMH1_dataMass0.0eV_fitMH1_fitMass0.0eV_dist10.0kpc_fit1D.txt")
    deltaT2D, nsig2D, nll2D = read_nll2T("fineNll2T_dataMH1_dataMass0.0eV_fitMH1_fitMass0.0eV_dist10.0kpc_fit2D.txt")
    plt.plot(deltaT1D, nll1D-np.min(nll1D), "o-")
    plt.plot(deltaT2D, nll2D-np.min(nll2D), "o-")
    plt.show()

