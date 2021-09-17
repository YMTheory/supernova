import numpy as np
import matplotlib.pyplot as plt
import ROOT
import uproot as up
import sys

ff = ROOT.TFile("out.root", "read")
hist1 = ROOT.Get("NO")
hist2 = ROOT.Get("IO")

nbinx, nbiny = hist1.GetNbinsX(), hist1.GetNbinsY()
lowx, higx = hist1.GetXaxis().GetBinCenter(1), hist1.GetXaxis().GetBinCenter(nbinx)
lowy, higy = hist1.GetYaxis().GetBinCenter(1), hist1.GetYaxis().GetBinCenter(nbiny)


cont = np.ones((nbiny, nbinx))
for i in range(nbinx):
    for j in range(nbiny):
        cont[j][i] = ( hist1.GetBinContent(i+1, j+1) - hist2.GetBinContent(i+1, j+1) ) 



fig = plt.figure(0)
im = plt.imshow(cont, origin="lower",interpolation="nearest",  vmin=0.01, extent=[lowx, higx, lowy, higy], aspect="auto")
fig.colorbar(im)
#plt.xlabel("time/s")
#plt.ylabel(r"$E_{dep}$/MeV")
#plt.ylabel(r"$E_{\nu}$/MeV")
#plt.savefig(sys.argv[4]+".pdf")

plt.show()
