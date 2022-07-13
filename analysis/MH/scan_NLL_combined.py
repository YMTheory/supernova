from tkinter import W
import numpy as np
import uproot as up
from iminuit import Minuit
import sys
import ROOT
import math
import matplotlib.pyplot as plt
plt.style.use("science")

#np.seterr(divide = 'ignore')
#np.seterr(divide = 'warn') 

chaname = ["eES", "pES", "IBD"]

def getIntegral(hist1d, x1, x2, opt):
    binx1 = hist1d.GetXaxis().FindBin(x1)
    binx2 = hist1d.GetXaxis().FindBin(x2)
    return hist1d.Integral(binx1, binx2, opt)


def calc_NLL(data, pdf, dT, Tstart, Tend)->float:
    nll = 0
    for i in data:
        mid = pdf.Interpolate(i+dT)
        if mid <= 0:
            mid = 1e-10
        nll += np.log(mid)
    nu = getIntegral(pdf, Tstart+dT, Tend+dT, "width")  # expected event numebr
    nll -= nu
    #N = len(data)

    #nll += N * np.log(nu)
    #nll -= nu

    return -nll



if __name__ == "__main__" :

    modName = None
    mod = None
    MH = "IO"
    ch = []
    Ethr = 0.2
    draw_flag = True

    if len(sys.argv) > 1:
        ### parametere configuration
        #print("argument number: %d"%int((len(sys.argv)-1)/2))
        for i in range(1, len(sys.argv)-1, 2):
            #print("==> Fitting ")
            if sys.argv[i] == "-modelName":
                modName = (sys.argv[i+1])
                #print("=====> ModelName: %s"%modName)
            elif sys.argv[i] == "-modelNo":
                mod = int(sys.argv[i+1])
                #print("=====> ModelNo: %d"%mod)
            elif sys.argv[i] == "-cha1" :
                ch.append(sys.argv[i+1])
                #print("=====> Channel 1: %s"%ch[-1])
            elif sys.argv[i] == "-cha2" :
                ch.append(sys.argv[i+1])
                #print("=====> Channel 2: %s"%ch[-1])
            elif sys.argv[i] == "-cha3" :
                ch.append(sys.argv[i+1])
                #print("=====> Channel 3: %s"%ch[-1])
            elif sys.argv[i] == "-mo" :
                MH = sys.argv[i+1]
                #print("=====> Mass ordering: %s"%MH)
            elif sys.argv[i] == "-Ethr" :
                Ethr = float(sys.argv[i+1])
                #print("=====> Ethr: %.2f MeV"%Ethr)
            elif sys.argv[i] == "-draw":
                draw_flag = sys.argv[i+1]
            else:
                print("Error: No such argument !")
                #exit(-1)


    outfn = "/junofs/users/miaoyu/supernova/analysis/MH/results/res_%s_%s_"%(modName+str(mod), MH)
    for cha in ch:
        outfn += cha
    outfn += "_%.2fMeV.txt"%Ethr
    with open(outfn,"w") as f:
        modArr = []
        if modName == "Garching" and mod == 32:   ## scanning all Garching models
            modArr = [81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502,92503,92700,92701,92702,92703,94000,94001,94002,94003]
        if modName == "Burrows2D" and mod == 6000:  ## scanning all Burrows 2D models
            modArr = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26]
        else:
            modArr = [mod]
        for mod in modArr:
            bestT1, bestT2 = [], []
            locMinNll1, locMinNll2 = [], []
    
            pdff1 = [None for i in range(len(ch))]
            pdff2 = [None for i in range(len(ch))]
            pdfh1 = [None for i in range(len(ch))]
            pdfh2 = [None for i in range(len(ch))]
            dataf = [None for i in range(len(ch))]
            pdf1_arr, pdf2_arr = [], []
            data_arr = []
            Nnu_arr = []
    
            for i, cha in enumerate(ch):
                ## read pdf
                if cha == "pES":
                    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/%s_PDF_NO_10kpc_%s_%.2fMeV.root" %(modName+str(mod), cha, Ethr)
                else:
                    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/%s_PDF_NO_10kpc_%s_%.2fMeV.root" %(modName+str(mod), cha, 0.20)
                hn = cha
                pdff1[i] = ROOT.TFile(pdffile, "read")
                #pdfh1[i] = pdff1[i].Get(hn)
                pdfh1[i] = pdff1[i].Get("h1")
                pdf1_arr.append(pdfh1[i])
                if cha == "pES":
                    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/%s_PDF_IO_10kpc_%s_%.2fMeV.root" %(modName+str(mod), cha, Ethr)
                else:
                    pdffile = "/junofs/users/miaoyu/supernova/production/PDFs/10kpc/%s_PDF_IO_10kpc_%s_%.2fMeV.root" %(modName+str(mod), cha, 0.20)
                hn = cha
                pdff2[i] = ROOT.TFile(pdffile, "read")
                pdfh2[i] = pdff2[i].Get("h1")
                pdf2_arr.append(pdfh2[i])
    
                # declare variables
                if cha == "pES":
                    datafile = "/junofs/users/miaoyu/supernova/production/Data/10kpc/%s_%s_data_%s_10kpc_thr%.2fMeV.root" %( modName+str(mod), cha, MH, Ethr)
                else:
                    datafile = "/junofs/users/miaoyu/supernova/production/Data/10kpc/%s_%s_data_%s_10kpc_thr%.2fMeV.root" %( modName+str(mod), cha, MH, 0.20)
                #print(datafile)
                ff = up.open(datafile)
                tt = ff["tFixedStat"]
                evtID = tt["evtID"].array()
                nuTime = tt["nuTime1D"].array()
                N_nu_per_sub = int(len(evtID)/500)
                #print("Neutrino events in each sub for %s %s channel with %.2f MeV threshold: %d" %(MH, ch, Ethr, N_nu_per_sub))
                Nnu_arr.append(N_nu_per_sub)
                data_arr.append(nuTime)
            
    
            ## time variable definitions : 
            Tmin, Tmax = 0, 60     # Fitting time range , unit : ms
            extTmin, extTmax = 10, 50
            dtmin, dtmax = -9, 9
            dT_arr = np.arange(dtmin, dtmax, 1)
            nstep = len(dT_arr)
                
            if draw_flag:
                ### Scanning
                fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(16, 4))
                stitle = "%s-%s-10kpc-%.2fMeV"%(modName+str(mod), MH, Ethr)
                for cha in ch:
                    stitle += '-'
                    stitle += cha
                st = fig.suptitle(stitle, fontsize="x-large")
    
            nStat = 500
            bestT1_arr, bestT2_arr, locMinNll1_arr, locMinNll2_arr = [], [], [], []
            for iStat in range(nStat):
                #if nStat % 10 == 0:
                #    print("Running SN event %d" %iStat)
                nll1, nll2 = np.zeros(nstep), np.zeros(nstep)
                for ich in range(len(ch)):
                    time_evt = data_arr[ich][iStat*Nnu_arr[ich]:(iStat+1)*Nnu_arr[ich]]
                    for j, dT in enumerate(dT_arr):
                        if modName == 'Garching':
                            Tstart, Tend = -20 , 20
                        if modName == "Burrows2D":
                            Tstart, Tend = 10, 50
                        tmp_nll1 = calc_NLL(time_evt, pdf1_arr[ich], dT, Tstart, Tend)
                        tmp_nll2 = calc_NLL(time_evt, pdf2_arr[ich], dT, Tstart, Tend)
                        nll1[j] += tmp_nll1
                        nll2[j] += tmp_nll2
    
                if draw_flag:
                    ax1.plot(dT_arr, 2*nll1, lw=2)    ## NO scanning
                    ax2.plot(dT_arr, 2*nll2, lw=2)    ## IO scanning
    
                bestT1, bestT2, locMinNll1, locMinNll2 = 0, 0, 10000, 10000
                for k, val in enumerate(nll1):
                    if val < locMinNll1:
                        locMinNll1 = val
                        bestT1 = dT_arr[k]
                for k, val in enumerate(nll2):
                    if val < locMinNll2:
                        locMinNll2 = val
                        bestT2 = dT_arr[k]
    
                locMinNll1_arr.append(locMinNll1)
                bestT1_arr.append(bestT1)
                locMinNll2_arr.append(locMinNll2)
                bestT2_arr.append(bestT2)
    
            locMinNll1_arr = np.array(locMinNll1_arr)
            bestT1_arr = np.array(bestT1_arr)
            locMinNll2_arr = np.array(locMinNll2_arr)
            bestT2_arr = np.array(bestT2_arr)
    
            if MH == "NO":
                f.write("%d mediumChi2 %.3f "%(mod, np.median(locMinNll2_arr - locMinNll1_arr)) + " ")
                for ic in range(len(Nnu_arr)):
                    f.write(str(Nnu_arr[ic]) + " ")
                f.write("\n")
            elif MH == "IO":
                f.write("%d mediumChi2 %.3f "%(mod, np.median(locMinNll1_arr - locMinNll2_arr)))
                for ic in range(len(Nnu_arr)):
                    f.write(str(Nnu_arr[ic]) + " ")
                f.write("\n")
    
    
            if draw_flag:
                ax1.set_xlabel(r"$\Delta T$ [ms]", fontsize=14)
                ax1.set_ylabel("-2ln(L)", fontsize=14)
                ax2.set_xlabel(r"$\Delta T$ [ms]", fontsize=14)
                ax2.set_ylabel("-2ln(L)", fontsize=14)
    
                ax3.hist(bestT1_arr, bins=20, range=(-10, 10), histtype="step", lw=2, label="NO fit: mean = %.2f ms"%np.mean(bestT1_arr))
                ax3.hist(bestT2_arr, bins=20, range=(-10, 10), histtype="step", lw=2, label="IO fit: mean = %.2f ms"%np.mean(bestT2_arr))
                ax3.set_xlabel(r"Best $\Delta T$ [ms]", fontsize=14)
                ax3.legend(prop={"size":14})
    
                if MH == "IO":
                    ax4.hist(2*locMinNll1_arr - 2*locMinNll2_arr, color="blue", edgecolor="black", bins=50, lw=2, label="medium = %.2f "%np.median(locMinNll1_arr - locMinNll2_arr))
                elif MH == "NO":
                    ax4.hist(2*locMinNll2_arr - 2*locMinNll1_arr, color="blue", edgecolor="black", bins=50, lw=2, label="medium = %.2f "%np.median(locMinNll2_arr - locMinNll1_arr))
                ax4.set_xlabel(r"$\Delta \chi^2$", fontsize=14)
                ax4.legend(prop={"size":14})
    
    
                plt.tight_layout()
                plt.savefig("/junofs/users/miaoyu/supernova/analysis/MH//results/"+stitle+".pdf")
                plt.show()
    



