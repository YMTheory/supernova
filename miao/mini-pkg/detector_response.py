import matplotlib.pyplot as plt
import numpy as np
import ROOT

libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)

imode = 82503
dist = 10

# === Configuration === #
snDet = ROOT.SNdetect.instance()
snDet.setSrcModel(imode)
modelSrc = snDet.getPointerSrc()

snDet.initChannel(1)
modelSrc.setSNDistance(dist)
Ethr = 0.1
snDet.getPointerEffectLS().setThresholdE(Ethr)
snDet.initFCN()

# === === === === === === #

def energy_spectra(MH):
    EArr = np.arange(0.1, 40, 0.1)
    T, Evis, Eobs = [], [], []
    for i in EArr:
        T.append(snDet.getTSpectrum(i, -1, MH))
        Evis.append(snDet.getEvisSpectrum(i, -1, MH))
        Eobs.append(snDet.getEobsSpectrum(i, -1, MH))
        print(i, T[-1], Evis[-1], Eobs[-1])

    EArr = np.array(EArr)
    T = np.array(T)
    Evis = np.array(Evis)
    Eobs = np.array(Eobs)

    return EArr, T, Evis, Eobs

def time_spectra(MH):
    timeArr = np.arange(-0.1, 0.1, 0.005)
    T, Evis, Eobs = [], [], []
    for i in timeArr:
        Evis.append(snDet.getEventAboveEthrVisAtTime(i, Ethr, -1, MH))
        Eobs.append(snDet.getEventAboveEthrObsAtTime(i, Ethr, -1, MH))
        print(i, Evis[-1], Eobs[-1])

    return timeArr, Evis, Eobs




if __name__ == "__main__":

    """
    timeArr, Evis, Eobs = time_spectra(1)
    plt.plot(timeArr, Evis, "-", label="w/o resolution")
    plt.plot(timeArr, Eobs, "-", label="w/  resolution")

    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("Event Rate")
    """

    EArr, T, Evis, Eobs = energy_spectra(1)
    plt.plot(EArr, T*EArr, "-", label="Edep")
    plt.plot(EArr, Evis*EArr, "-", label="w/o resolution")
    plt.plot(EArr, Eobs*EArr, "-", label="w/  resolution")

    plt.legend()
    plt.semilogx()
    plt.semilogy()
    plt.xlabel("energy/MeV")
    plt.ylabel("Event Rate")
    plt.show()




