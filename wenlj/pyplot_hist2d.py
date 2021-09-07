import numpy as np
import matplotlib.pyplot as plt
import ROOT
import uproot as up
import sys
import boost_histogram as bh




if sys.argv[1] == "single":


    if len(sys.argv) != 5:
        print("Wrong parameter number inputs!")
        sys.exit(1)
    
    fin = ROOT.TFile(sys.argv[2], "read")
    hist = fin.Get(sys.argv[3])    # hist2d
    nbinx, nbiny = hist.GetNbinsX(), hist.GetNbinsY()
    lowx, higx = hist.GetXaxis().GetBinCenter(1), hist.GetXaxis().GetBinCenter(nbinx)
    lowy, higy = hist.GetYaxis().GetBinCenter(1), hist.GetYaxis().GetBinCenter(nbiny)
    wx, wy = hist.GetXaxis().GetBinWidth(0), hist.GetYaxis().GetBinWidth(0)
    
    # projection on X, Y
    hTime = hist.ProjectionX()
    hE    = hist.ProjectionY()
    cent_time, cent_E = [], []
    cont_time, cont_E = [], []
    width_time, width_E = [], []
    for i in range(nbinx):
        cent_time.append(hTime.GetBinCenter(i))
        cont_time.append(hTime.GetBinContent(i))
        width_time.append(hTime.GetBinWidth(i))
    
    for j in range(nbiny):
        cent_E.append(hE.GetBinCenter(j))
        cont_E.append(hE.GetBinContent(j))
        width_E.append(hE.GetBinWidth(i))
    
    cent_time  = np.array(cent_time)
    cont_time  = np.array(cont_time)
    width_time = np.array(width_time)
    cent_E  = np.array(cent_E)
    cont_E  = np.array(cont_E)
    width_E = np.array(width_E)
    
    
    #xbins, ybins =[], []
    cont = np.ones((nbiny, nbinx))
    for i in range(nbinx):
        for j in range(nbiny):
            cont[j][i] = hist.GetBinContent(i+1, j+1)
            #xbins.append(hist.GetXaxis().GetBinCenter(i+1)-wx)
            #ybins.append(hist.GetYaxis().GetBinCenter(i+1)-wy)
    
    #xbins.append(hist.GetXaxis().GetBinCenter(nbinx)+wx)
    #ybins.append(hist.GetYaxis().GetBinCenter(nbiny)+wy)
    
    
    
    # 2D energy-time spectrum
    fig = plt.figure(0)
    im = plt.imshow(cont, origin="lower",interpolation="nearest",  vmin=0.01, extent=[lowx, higx, lowy, higy], aspect="auto")
    fig.colorbar(im)
    plt.xlabel("time/s")
    #plt.ylabel(r"$E_{dep}$/MeV")
    plt.ylabel(r"$E_{\nu}$/MeV")
    plt.savefig(sys.argv[4]+".pdf")
    
    plt.figure(1)
    plt.plot(cent_time, cont_time*width_E[0], "-")
    plt.xlabel("time/ns")
    plt.ylabel("Event Rate")
    plt.savefig(sys.argv[4]+"_ProjT.pdf")
    
    
    plt.figure(2)
    plt.plot(cent_E, cont_E*width_time[0], "-")
    plt.semilogx()
    plt.semilogy()
    plt.xlim(0.2, 10)
    plt.xlabel("E/MeV")
    plt.ylabel("Event Rate")
    plt.savefig(sys.argv[4]+"_ProjE.pdf")
    plt.show()





if sys.argv[1] == "multi":
    totFile = int(sys.argv[2])

    Tmin, Tmax = -0.04, 0.1
    Emin, Emax = 0.1, 25
    Tmin, Tmax = -0.1, 0.1
    Emin, Emax = 0.1, 10

    filelist, histlist, outlist = [], [], []
    for i in range(totFile):
        filelist.append(sys.argv[3+i])
        histlist.append(sys.argv[3+totFile+i])
    for i in range(1):
        outlist.append(sys.argv[3+2*totFile+i])


    cent_time, cent_E = [[] for i in range(totFile)], [[] for i in range(totFile)]
    cont_time, cont_E = [[] for i in range(totFile)], [[] for i in range(totFile)]
    width_time, width_E = [[] for i in range(totFile)], [[] for i in range(totFile)]

    for i in range(totFile):
        fin = ROOT.TFile(filelist[i], "read")
        hist = fin.Get(histlist[i])    # hist2d
        nbinx, nbiny = hist.GetNbinsX(), hist.GetNbinsY()
        lowx, higx = hist.GetXaxis().GetBinCenter(1), hist.GetXaxis().GetBinCenter(nbinx)
        lowy, higy = hist.GetYaxis().GetBinCenter(1), hist.GetYaxis().GetBinCenter(nbiny)

        binx1 = hist.GetXaxis().FindBin(Tmin)
        binx2 = hist.GetXaxis().FindBin(Tmax)
        biny1 = hist.GetYaxis().FindBin(Emin)
        biny2 = hist.GetYaxis().FindBin(Emax)
        print(binx1, binx2, biny1, biny2, hist.Integral(binx1, binx2, biny1, biny2, "width"))
        

        # projection on X, Y
        hTime = hist.ProjectionX()
        hE    = hist.ProjectionY()
        for j in range(nbinx):
            cent_time[i].append(hTime.GetBinCenter(j+1))
            cont_time[i].append(hTime.GetBinContent(j+1))
            width_time[i].append(hE.GetBinWidth(j+1))
        
        for j in range(nbiny):
            cent_E[i].append(hE.GetBinCenter(j+1))
            cont_E[i].append(hE.GetBinContent(j+1))
            width_E[i].append(hTime.GetBinWidth(j+1))
        
    
    cent_time  = np.array(cent_time)
    cont_time  = np.array(cont_time)
    width_time = np.array(width_time)
    cent_E  = np.array(cent_E)
    cont_E  = np.array(cont_E)
    width_E = np.array(width_E)

    #labels = ["no osc", "NO", "IO"] 
    labels = ["0.0eV", "0.6eV", "1.0eV"]

    plt.figure(0)
    for i in range(totFile):
        print(np.sum(cont_time[i]*width_time[i]))
        plt.plot(cent_time[i], cont_time[i]*width_time[i]/50, "-", label=labels[i])
    plt.xlabel("time/s")
    plt.ylabel(r"dN/dt ($s^{-1}$)")
    plt.legend()
    plt.savefig(outlist[0]+"_ProjT.pdf")
    
    
    plt.figure(1)
    for i in range(totFile):
        plt.plot(cent_E[i], cont_E[i]*width_E[i]*cent_E[i], "-", label=labels[i])
    plt.xlabel("E/MeV")
    plt.ylabel("E dN/dE")
    plt.legend()
    plt.semilogx()
    plt.semilogy()
    plt.savefig(outlist[0]+"_ProjE.pdf")
    
    #plt.show()
        








