import numpy as np
import uproot as up
from iminuit import Minuit
import sys
import ROOT
import math
import matplotlib.pyplot as plt
plt.style.use("science")

do_fine_scan = False
combined_pES_time = True

def getIntegral(hist1d, x1, x2, opt):
    binx1 = hist1d.GetXaxis().FindBin(x1)
    binx2 = hist1d.GetXaxis().FindBin(x2)
    return hist1d.Integral(binx1, binx2, opt)

def getIntegral(hist1d, x1, x2, opt):
    binx1 = hist1d.GetXaxis().FindBin(x1)
    binx2 = hist1d.GetXaxis().FindBin(x2)
    return hist1d.Integral(binx1, binx2, opt)


def calc_NLL(data, pdf, dT)->float:
    nll = 0
    for i in data:
        nll += np.log(pdf.Interpolate(i+dT))
    nu = getIntegral(pdf, -20+dT, 20+dT, "width")  # expected event numebr
    nll -= nu
    N = len(data)

    nll += N * np.log(nu)
    nll -= nu

    return -nll



if __name__ == "__main__" :

    mod = 82503
    MH = "IO"
    cha = "pES"
    Ethr = 0.2

    if len(sys.argv) > 1:
        ### parametere configuration
        print("argument number: %d"%int((len(sys.argv)-1)/2))
        for i in range(1, len(sys.argv)-1, 2):
            print("==> Fitting ")
            if sys.argv[i] == "-model":
                mod = int(sys.argv[i+1])
                print("=====> Model: %d"%mod)
            elif sys.argv[i] == "-cha" :
                cha = sys.argv[i+1]
                print("=====> Channel: %s"%cha)
            elif sys.argv[i] == "-mo" :
                MH = sys.argv[i+1]
                print("=====> Mass ordering: %s"%MH)
            elif sys.argv[i] == "-Ethr" :
                Ethr = float(sys.argv[i+1])
                print("=====> Ethr: %.2f MeV"%Ethr)
            else:
                print("Error: No such argument !")
                exit(-1)



    if combined_pES_time:
        fig0, ax0 = plt.subplots(figsize=(6, 6))
        time0_pES = []

        pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Garching%d_PDF_NO_10kpc_%s_%.2fMeV.root" %(mod, "pES", 0.15)
        hn = "pES"
        pdff = ROOT.TFile(pdffile, "read")
        pdf = pdff.Get(hn)

        Tmin, Tmax = -20, 20     # Fitting time range , unit : ms
        extTmin, extTmax = -30, 30
        datafile = "/junofs/users/miaoyu/supernova/production/Data/10kpc/Garching%d_%s_data_%s_10kpc_thr%.2fMeV.root" %(mod, "pES", "IO", Ethr)
        print(datafile)
        ff = up.open(datafile)
        tt = ff["tFixedStat"]
        evtID = tt["evtID"].array()
        nuTime = tt["nuTime"].array()
        N_nu_per_sub = int(len(evtID)/500)
        print("Neutrino events in each sub: %d" %N_nu_per_sub)
        dtmin, dtmax = -9, 10
        dT_arr = np.arange(dtmin, dtmax, 1)


        nStat = 500
        for iSub in range(nStat):
            if iSub % 10 == 0:
                print("Running sub %d" %iSub)
            # calc NLL:
            time_evt = nuTime[iSub*N_nu_per_sub:(iSub+1)*N_nu_per_sub]
            nll_arr = []
            for dT in dT_arr:
                nll = calc_NLL(time_evt, pdf, dT)
                nll_arr.append(nll)
        
            locNLLMin, locT = 100, 0
            for i in range(len(nll_arr)):
                if nll_arr[i] < locNLLMin:
                    locNLLMin = nll_arr[i]
                    locT1 = dT_arr[i]

            time0_pES.append(locT1)

        ax0.hist(time0_pES, bins=20, range=(-10, 10), color="blue", edgecolor="black")
        ax0.set_xlabel(r"$t_0$ [ms]", fontsize=14)


    ## read pdf
    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Garching%d_PDF_NO_10kpc_%s_%.2fMeV.root" %(mod, cha, Ethr)
    hn = cha
    pdff1 = ROOT.TFile(pdffile, "read")
    pdf1 = pdff1.Get(hn)
    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Garching%d_PDF_IO_10kpc_%s_%.2fMeV.root" %(mod, cha, Ethr)
    hn = cha
    pdff2 = ROOT.TFile(pdffile, "read")
    pdf2 = pdff2.Get(hn)

    # declare variables
    Tmin, Tmax = -20, 20     # Fitting time range , unit : ms
    extTmin, extTmax = -30, 30
    datafile = "/junofs/users/miaoyu/supernova/production/Data/10kpc/Garching%d_%s_data_%s_10kpc_thr%.2fMeV.root" %(mod, cha, MH, Ethr)
    print(datafile)
    ff = up.open(datafile)
    tt = ff["tFixedStat"]
    evtID = tt["evtID"].array()
    nuTime = tt["nuTime"].array()
    N_nu_per_sub = int(len(evtID)/500)
    print("Neutrino events in each sub: %d" %N_nu_per_sub)
    dtmin, dtmax = -9, 10
    dT_arr = np.arange(dtmin, dtmax, 1)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(16, 4))
    nStat = 500
    nllMin1, nllMin2 = [], []
    bestT1, bestT2 = [], []
    for iSub in range(nStat):
        if iSub % 10 == 0:
            print("Running sub %d" %iSub)
        # calc NLL:
        time_evt = nuTime[iSub*N_nu_per_sub:(iSub+1)*N_nu_per_sub]
        nll_arr1, nll_arr2 = [], []
        nll_arr_fine1, nll_arr_fine2 = [], []
        tarr_fine1, tarr_fine2 = [], []
        for dT in dT_arr:
            nll1 = calc_NLL(time_evt, pdf1, dT)
            nll2 = calc_NLL(time_evt, pdf2, dT)
            nll_arr1.append(nll1)
            nll_arr2.append(nll2)
    
        if not do_fine_scan:
            ax1.plot(dT_arr, nll_arr1, lw=1.5)
            ax2.plot(dT_arr, nll_arr2, lw=1.5)

        locNLLMin1, locNLLMin2, locT1, locT2 = 100, 100, 0, 0
        for i in range(len(nll_arr1)):
            if nll_arr1[i] < locNLLMin1:
                locNLLMin1 = nll_arr1[i]
                locT1 = dT_arr[i]
            if nll_arr2[i] < locNLLMin2:
                locNLLMin2 = nll_arr2[i]
                locT2 = dT_arr[i]

        if do_fine_scan:
            tstep = 0.1
            for istep in range(21):
                ddT = tstep * (istep-10)
                if locT1 + ddT > extTmax or locT1 + ddT < extTmin or locT2 + ddT > extTmax or locT2 + ddT < extTmin:
                    print(locT1, locT1+ddT, locT2, locT2+ddT)
                    print("Warning! time is shifted outside the pdf!!!")
                    continue
                nll = calc_NLL(time_evt, pdf1, ddT+locT1)
                nll_arr_fine1.append(nll)
                tarr_fine1.append(ddT+locT1)

                nll = calc_NLL(time_evt, pdf2, ddT+locT2)
                nll_arr_fine2.append(nll)
                tarr_fine2.append(ddT+locT2)
                
            locNLLMin1, locNLLMin2, locT1, locT2 = 100, 100, 0, 0
            for i in range(len(nll_arr_fine1)):
                if nll_arr_fine1[i] < locNLLMin1:
                    locNLLMin1 = nll_arr_fine1[i]
                    locT1 = tarr_fine1[i]
                if nll_arr_fine2[i] < locNLLMin2:
                    locNLLMin2 = nll_arr_fine2[i]
                    locT2 = tarr_fine2[i]
        
            ax1.plot(tarr_fine1, nll_arr_fine1, lw=1.5)
            ax2.plot(tarr_fine2, nll_arr_fine2, lw=1.5)

        bestT1.append(locT1)
        nllMin1.append(locNLLMin1)
        bestT2.append(locT2)
        nllMin2.append(locNLLMin2)

    nllMin1 = np.array(nllMin1)
    nllMin2 = np.array(nllMin2)
    bestT1 = np.array(bestT1)
    bestT2 = np.array(bestT2)
    
    ax3.hist(bestT1, bins=20, range=(-10, 10), histtype="step", label="NO fit : %.2f"%np.mean(bestT1))
    ax3.hist(bestT2, bins=20, range=(-10, 10), histtype="step", label="NO fit : %.2f"%np.mean(bestT2))

    #ax4.hist(nllMin1, bins=50, histtype="step", label="NO fit : %.2f"%np.mean(nllMin1))
    #ax4.hist(nllMin2, bins=50, histtype="step", label="IO fit : %.2f"%np.mean(nllMin2))
    if MH == "NO":
        ax4.hist(nllMin2 - nllMin1, bins=50, label="mean = %.2f"%np.mean(nllMin2 - nllMin1))
    elif MH == "IO":
        ax4.hist(nllMin1 - nllMin2, bins=50, label="mean = %.2f"%np.mean(nllMin1 - nllMin2))

    ax1.set_xlabel(r"$\Delta T$ [ms]", fontsize=14)
    ax1.set_ylabel("NLL", fontsize=14)
    ax2.set_xlabel(r"$\Delta T$ [ms]", fontsize=14)
    ax2.set_ylabel("NLL", fontsize=14)
    ax3.set_xlabel("best dT [ms]", fontsize=14)
    ax3.legend(loc="upper right", prop={'size':14})
    ax4.set_xlabel(r"$\Delta$NLL", fontsize=14)
    ax4.legend(loc="upper right", prop={'size':14})

    plt.tight_layout()
    plt.show()






