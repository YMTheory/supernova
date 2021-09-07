import numpy as np
import matplotlib.pyplot as plt
import ROOT

def loadTH2D(filename, histname):
    ff = ROOT.TFile(filename, "read")
    hh = ff.Get(histname)
    xbins, ybins = [], []
    cont = np.zeros((hh.GetNbinsY(), hh.GetNbinsX()))
    
    for i in range(hh.GetNbinsX()):
        for j in range(hh.GetNbinsY()):
            cont[j, i] = hh.GetBinContent(i+1, j+1)

    return cont
        




import namespace as ns

def merge():

    mod = 82503
    dist = 2
    #nuMass = 0.1
    nuType = -1
    chaname = 1
    MH = 1
    
    # NuE channel
    Evismin = np.arange(0.0, 25.0, 0.5)
    Evismax = np.arange(0.5, 25.5, 0.5)
    Evismin[0] = 0.1
    # NuP channel
    #Evismin = np.arange(0.1, 5.0, 0.1)
    #Evismax = np.arange(0.2, 5.1, 0.1)
    #Evismin[0] = 0.1
    Emin, Emax = Evismin[0], Evismax[-1]


    plotMass = 12

    for num in range(0, 20, 1):
        
        cont0, cont1, cont2 = [], [], []

        nuMass = num/10.
        print("Current nuMass : %.2f"%nuMass)
        for no in range(len(Evismin)):
            emin, emax = Evismin[no], Evismax[no]
            rootfile = ns.pdfFileName(mod, dist, chaname, nuMass, MH, nuType, emax, no)
            
            hh0 = ns.pdfEvisHist2DName(mod, chaname, MH)
            cont0_tmp = loadTH2D(rootfile, hh0)
            cont0.append(cont0_tmp)

            hh1 = ns.pdfEvHist2DName(mod, chaname, MH)
            cont1_tmp = loadTH2D(rootfile, hh1)
            cont1.append(cont1_tmp)


            hh2 = ns.pdfBkgHist2DName(mod, chaname, MH)
            cont2_tmp = loadTH2D(rootfile, hh2)
            cont2.append(cont2_tmp)
            
        res0 = np.concatenate(cont0, axis=0)
        res1 = np.zeros(cont1[0].shape)
        for arr in cont1:
            res1 = res1 + arr
        res2 = np.concatenate(cont2, axis=0)



        out = ROOT.TFile(ns.pdfSumFileName(mod, dist, chaname, nuMass, MH, nuType, emax), "recreate")
        hist0 = ROOT.TH2D("TEvisPdf" , "TEvisPdf", res0.shape[1], -0.1, 0.1, res0.shape[0], Emin, Emax)
        hist1 = ROOT.TH2D("TEnuPdf", "TEnuPdf", res1.shape[1], -0.1, 0.1, res1.shape[0], 0, 60)
        hist2 = ROOT.TH2D("TBkgEvisPdf", "TBkgEvisPdf", res2.shape[1], -0.1, 0.1, res2.shape[0], Emin, Emax)


        for i in range(res0.shape[1]):
            for j in range(res0.shape[0]):
                hist0.SetBinContent(i+1, j+1, res0[j, i])

        for i in range(res1.shape[1]):
            for j in range(res1.shape[0]):
                hist1.SetBinContent(i+1, j+1, res1[j, i])

        for i in range(res2.shape[1]):
            for j in range(res2.shape[0]):
                hist2.SetBinContent(i+1, j+1, res2[j, i])



        hist0.Write()
        hist1.Write()
        hist2.Write()

        out.Close()

        if nuMass == plotMass :

            fig = plt.figure(0)
            im = plt.imshow(res0, origin="lower",interpolation="nearest",  vmin=0.01, extent=[-0.1, 0.1, 0.1, 25], aspect="auto")
            fig.colorbar(im)
            plt.xlabel("time/s")
            plt.ylabel(r"$E_{dep}$/MeV")
            plt.title("signal spectra")
             
            fig = plt.figure(1)
            im = plt.imshow(res1, origin="lower",interpolation="nearest",  vmin=0.01, extent=[-0.1, 0.1, 0.0, 60], aspect="auto")
            fig.colorbar(im)
            plt.xlabel("time/s")
            plt.ylabel(r"$E_{\nu}$/MeV")
            plt.title("neutrino spectra")
             
            fig = plt.figure(2)
            im = plt.imshow(res2, origin="lower",interpolation="nearest",  vmin=0.01, extent=[-0.1, 0.1, 0.1, 25], aspect="auto")
            fig.colorbar(im)
            plt.xlabel("time/s")
            plt.ylabel(r"$E_{dep}$/MeV")
            plt.title("background spectra")
             

            plt.show()


if __name__ == "__main__":
    merge()



