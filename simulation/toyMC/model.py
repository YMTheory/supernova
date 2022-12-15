import numpy as np
import matplotlib.pyplot as plt
import uproot as up
from scipy.special import gamma


def model(a, T, t0, t):
    if t0 < t < t0 + T:
        return 1 + a * np.sin(t)
    else:
        return 1



def generator(a, T, t0):

    t = np.arange(0, 1000, 0.1)
    Nevt = 10000
    Npoint = 1000
    y = np.zeros((Nevt, Npoint))

    for ievt in range(Nevt):
        for ipt in range(Npoint):
            t = np.random.uniform
            y[ievt, ipt] = model(a, T, t0, t)
    


def toyest_model(muNO, muIO):
    sens_asimov_NO = -2 * (muNO * np.log(muIO) - muIO - np.log(gamma(muNO+1)))
    nNO = np.random.poisson(muNO, size=10000)
    sens_toyMC_NO = []
    for i in nNO:
        if i != 0:
            tmp_sens =  -2 * (i * np.log(muIO) - muIO - np.log(gamma(i+1)))
        else:
            tmp_sens = 2 * muIO

        sens_toyMC_NO.append(tmp_sens)
    sens_toyMC_NO = np.array(sens_toyMC_NO)

    return sens_asimov_NO, sens_toyMC_NO





if __name__ == "__main__":

    muNO = np.random.uniform(0, 1, size=100)
    muIO = np.random.uniform(0, 1, size=100)
    
    sens_asimov_NO = 0
    sens_toyMC_NO = np.zeros(10000)
    for i, j in zip(muNO, muIO):
        tmp_sens_asimov_NO, tmp_sens_toyMC_NO = toyest_model(i, j)
        sens_toyMC_NO += tmp_sens_toyMC_NO
        sens_asimov_NO += tmp_sens_asimov_NO

    print(f"Asimov data sensitivity = {sens_asimov_NO}")
    print(f"ToyMC data median sensitivity = {np.median(sens_toyMC_NO)}")
    print(f"ToyMC data mean sensitivity = {np.mean(sens_toyMC_NO)}")
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(np.arange(0, 100, 1), muNO, label="NO")
    ax1.plot(np.arange(0, 100, 1), muIO, label="IO")
    ax1.legend()
    ax0.hist(sens_toyMC_NO, bins=50, histtype="step")
    ax0.vlines(sens_asimov_NO, 0, 600, color="red", lw=2, label="Asimov")
    ax0.vlines(np.median(sens_toyMC_NO), 0, 600, color="orange", lw=2, linestyles=":", label="ToyMC median")
    ax0.vlines(np.mean(sens_toyMC_NO), 0, 600, color="green", lw=2, linestyles="-.", label="ToyMC mean")
    ax0.legend()
    plt.show()






