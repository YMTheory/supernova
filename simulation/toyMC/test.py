from channel_analyser import channel
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

IBD = channel("IBD", "IO", "Garching", 82703, 0.15)
IBD.setNevtPerFile(100000)
IBD.setStartEvtId(2000)
IBD.setEndEvtId(2500)
#IBD.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/scale10/Garching82703_IBD_data_NO_10kpc_thr0.15MeV_Tmin10msTmax50ms_merger.root")
IBD._load_data()
IBD._load_pdf()


if 0:
    fig, ax = plt.subplots()
    t = np.arange(0, 50, 1)
    ax.plot(t, IBD._pdfNO_func(t), lw=2, label="NO")
    ax.plot(t, IBD._pdfIO_func(t), lw=2, label="IO")
    ax.set_xlabel("time [ms]", fontsize=15)
    ax.set_ylabel(r"d$N$/dt [ms$^{-1}$]", fontsize=15)
    ax.tick_params(axis="both", labelsize=14)
    ax.legend(prop={"size":14})
    plt.tight_layout()
    plt.show()


Nexp = 500
tmax, tmin = 50, 0
stepT = 0.01
Nstep = int((tmax - tmin) / stepT)

print(f"Total binning number = {Nstep}")

#######################################################################################################
############################ Asimov Dataset Sensitivity. ##############################################
nll_asimov = 0
for isig in range(Nstep):
    tt = tmin + isig * stepT
    n = IBD._pdfIO_func(tt) * stepT
    s = IBD._pdfNO_func(tt) * stepT 
    if s != 0:
        nll_asimov += -2 * (n * np.log(s) - s - np.log(gamma(n+1)))
    if s == 0:
        nll_asimov += 2 * s
    s = IBD._pdfIO_func(tt) * stepT
    if s != 0:
        nll_asimov += 2 * (n * np.log(s) - s - np.log(gamma(n+1)))
    if s == 0:
        nll_asimov -= 2 * s

print(f"Asimov dataset sensitivity = {nll_asimov:.2f}")

#######################################################################################################
################################# ToyData Sensitivity. ################################################
nll_arr = []
for iexp in range(Nexp):
    signals = np.zeros(Nstep) # generated toy data bin by bin
    for isig in range(Nstep):
        n = IBD._pdfIO_func(tmin + stepT * isig) * stepT
        s = np.random.poisson(n, size=1)
        if s <= 0:
            s = 0
        signals[isig] = s
    
    #arr = IBD.get_one_event(iexp+2000)
    #signals, _ = np.histogram(arr, bins=Nstep, range=(tmin, tmax))
    print(signals)
    print(f"===========Exp {iexp} Running, total {np.sum(signals):.2f} neutrinos===========")
    nll = 0
    for isig in range(Nstep):
        n = signals[isig]
        s = IBD._pdfNO_func(tmin + stepT * isig) * stepT
        if s != 0:
            nll += -2 * (n * np.log(s) - s - np.log(gamma(n+1)))
        if s == 0:
            nll += 2 * s

        s = IBD._pdfIO_func(tmin + stepT * isig) * stepT
        if s != 0:
            nll += 2 * (n * np.log(s) - s - np.log(gamma(n+1)))
        if s == 0:
            nll -= 2 * s
        
    nll_arr.append(nll)
#######################################################################################################
#######################################################################################################
    
    
fig, ax = plt.subplots()
ax.hist(nll_arr, bins=50, histtype="step")
ax.vlines(np.median(nll_arr), 0, 30, lw=2, linestyle=":", color="green", label=f"ToyMC median: {np.median(nll_arr):.2f}")
ax.vlines(np.mean(nll_arr), 0, 30, lw=2, linestyles="--", color="orange", label=f"ToyMC mean: {np.mean(nll_arr):.2f}")
ax.vlines(nll_asimov, 0, 30, lw=2, color="red", label=f"Asimov: {nll_asimov:.2f}")
ax.legend(prop={"size":14})
ax.set_xlabel(r"\Delta \chi^2", fontsize=14)
plt.show()
        
    
    



























