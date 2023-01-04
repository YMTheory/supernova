import numpy as np
import matplotlib.pyplot as plt
import ROOT
from channel_analyser import channel

def load_C14():
    f = ROOT.TFile("/junofs/users/miaoyu/supernova/production/PDFs/backgrounds/C14/C14_rate.root", "read")
    glow = f.Get("c14_low")
    ghig = f.Get("c14_high")

    return glow, ghig


def generate_C14(level, glow, ghig):
    if level == "low":
        c14rate = glow.Eval(Ethr)
    elif level == "high":
        c14rate = ghig.Eval(Ethr)
    else:
        print("Unknown level description!")

    nc14 = c14rate * (fitTmax - fitTmin)* 1e-3
    nc14 = np.random.poisson(nc14)
    if nc14 <= 0:
        return []
    Tc14 = np.random.uniform(fitTmin, fitTmax, size=nc14)
    print(f"C14-rate = {c14rate} Hz with event number in window = {nc14}.")
    return Tc14




if __name__ == "__main__" :

    MO      = "IO"
    model   = "Garching"
    modelNo = 82703
    Ethr    = 0.15
    fitTmin = -20
    fitTmax = 20
    fileNo  = 0
    dist    = 10
    exp     = "JUNO"
    level   = "high"
    startevt= 0
    endevt  = 10

    glow, ghig = load_C14()
    if level == "high":
        c14rate = ghig.Eval(Ethr)
    elif level == "low":
        c14rate = glow.Eval(Ethr)

    channels = {}
    channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    channels["pES"].c14rate = c14rate
    channels["IBD"] = channel("IBD", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)
    channels["eES"] = channel("eES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax=fitTmax, fileNo=fileNo, dist=dist, exp=exp)


    for cha in channels.values():
        cha.setNevtPerFile(1e5)
        cha.setStartEvtId(startevt)
        cha.setEndEvtId(endevt)
        cha.setDataFilePath(f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha.name}_binneddata_{MO}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_merger.root")
        cha._load_data_ak()
    






