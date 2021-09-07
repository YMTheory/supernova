import numpy as np
import matplotlib.pyplot as plt
import ROOT
import os, sys

def projectionX(filename):

    ff = ROOT.TFile(filename, "read")
    if not os.path.exists(filename):
        print("%s does not exist!"%filename)
        sys.exit(-1)
    inputHist2D = ff.Get("TEvisPdf")

    nbinsX = inputHist2D.GetNbinsX()
    nbinsY = inputHist2D.GetNbinsY()

    cent = np.zeros(nbinsX)
    cont = np.zeros(nbinsX)

    for i in range(1, nbinsX+1, 1):
        tmp_cont = 0.
        for j in range(1, nbinsY+1, 1):
            tmp_cont += inputHist2D.GetBinContent(i, j) * inputHist2D.GetYaxis().GetBinWidth(j)
        cent[i-1] = inputHist2D.GetXaxis().GetBinCenter(i)
        cont[i-1] = tmp_cont


    return cent, cont



def projectionY(filename):

    ff = ROOT.TFile(filename, "read")
    if not os.path.exists(filename):
        print("%s does not exist!"%filename)
        sys.exit(-1)
    inputHist2D = ff.Get("TEvisPdf")

    nbinsX = inputHist2D.GetNbinsX()
    nbinsY = inputHist2D.GetNbinsY()

    cent = np.zeros(nbinsY)
    cont = np.zeros(nbinsY)

    for j in range(1, nbinsY+1, 1):
        tmp_cont = 0.
        for i in range(1, nbinsX+1, 1):
            tmp_cont += inputHist2D.GetBinContent(i, j) * inputHist2D.GetXaxis().GetBinWidth(i)
        cent[j-1] = inputHist2D.GetYaxis().GetBinCenter(j)
        cont[j-1] = tmp_cont

    return cent, cont




if __name__ == "__main__":

    filename1 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/10kpc/NO/TEvisPDF_mod82503_cha1nue_nuType-1_mh1_mNu0.0eV_10.0kpc_0.1s_Evismax25.00_Sum.root"
    filename2 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/5kpc/NO/TEvisPDF_mod82503_cha1nue_nuType-1_mh1_mNu0.0eV_5.0kpc_0.1s_Evismax25.00_Sum.root"
    filename3 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/2kpc/NO/TEvisPDF_mod82503_cha1nue_nuType-1_mh1_mNu0.0eV_2.0kpc_0.1s_Evismax25.00_Sum.root"

    centTNO10kpc, contTNO10koc = projectionX(filename1)
    centTNO5kpc,  contTNO5kpc  = projectionX(filename2)
    centTNO2kpc,  contTNO2kpc  = projectionX(filename3)


    plt.plot(centTNO10kpc, contTNO10koc, "-",  label="NO 10kpc")
    plt.plot(centTNO5kpc, contTNO5kpc, "--",   label="NO 5kpc")
    plt.plot(centTNO2kpc, contTNO2kpc, "--",   label="NO 2kpc")

    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"Event Rate [$s^{-1}$]")
    plt.show()











