import numpy as np
import matplotlib.pyplot as plt
import ROOT

def Integral(filename):
    ff = ROOT.TFile(filename, "read")
    inputHist2D = ff.Get("TEvisPdf")

    Tmin, Tmax = -0.015, 0.02
    Emin, Emax = 0.1, 4

    x1 = inputHist2D.GetXaxis().FindBin(Tmin)
    x2 = inputHist2D.GetXaxis().FindBin(Tmax)
    y1 = inputHist2D.GetYaxis().FindBin(Emin)
    y2 = inputHist2D.GetYaxis().FindBin(Emax)

    return inputHist2D.Integral(x1, x2, y1, y2, "width")





if __name__ == "__main__" :

    filename1 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/10kpc/NO/TEvisPDF_mod82503_cha1nue_nuType-1_mh1_mNu0.0eV_10.0kpc_0.1s_Evismax25.00_Sum.root"
    filename2 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/10kpc/IO/TEvisPDF_mod82503_cha1nue_nuType-1_mh2_mNu0.0eV_10.0kpc_0.1s_Evismax25.00_Sum.root"

    #print(Integral(filename1))
    print(Integral(filename2))



    filename1 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/5kpc/NO/TEvisPDF_mod82503_cha1nue_nuType-1_mh1_mNu0.0eV_5.0kpc_0.1s_Evismax25.00_Sum.root"
    filename2 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/5kpc/IO/TEvisPDF_mod82503_cha1nue_nuType-1_mh2_mNu0.0eV_5.0kpc_0.1s_Evismax25.00_Sum.root"

    print(Integral(filename1))
    print(Integral(filename2))



    filename1 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/2kpc/NO/TEvisPDF_mod82503_cha1nue_nuType-1_mh1_mNu0.0eV_2.0kpc_0.1s_Evismax25.00_Sum.root"
    filename2 = "/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/2kpc/IO/TEvisPDF_mod82503_cha1nue_nuType-1_mh2_mNu0.0eV_2.0kpc_0.1s_Evismax25.00_Sum.root"

    print(Integral(filename1))
    print(Integral(filename2))















