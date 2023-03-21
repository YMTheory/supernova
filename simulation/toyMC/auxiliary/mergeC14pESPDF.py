import sys
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from array import array

out = True
if len(sys.argv) > 1:
    Ethr = float(sys.argv[1])

#pESfilename = "./Garching82703_nuePDF_NO_10kpc_pES_nuMass0.0eV_TEobs2dPDF_v2.root"
pESfilename = "./Garching82703_nuePDF_NO_10kpc_pES_nuMass0.0eV_TEobs2dPDF_JUNO_rebin.root"
c14filename = "/junofs/users/miaoyu/supernova/production/PDFs/backgrounds/C14/C14_visible_spectrum.root"

pESfile = ROOT.TFile(pESfilename, "read")
pESHist = pESfile.Get("h1")
Ntime = pESHist.GetNbinsX()
Nenergy = pESHist.GetNbinsY()
print(f"X-bin # {Ntime}, Y-bin # {Nenergy}.")

c14file = ROOT.TFile(c14filename, "read")
c14Hist = c14file.Get("c14")

## calculate C14 rate:
c14rate0 = 3.54e2   # 1e-18g/g, unit Hz -> 0.10 MeV
binmin0 = c14Hist.FindBin(0.10)
itot0 = c14Hist.Integral(binmin0, c14Hist.GetNbinsX())
binmin = c14Hist.FindBin(Ethr)
itot = c14Hist.Integral(binmin, c14Hist.GetNbinsX())
c14rate = itot / itot0 * c14rate0
nC14persec = c14rate * 1
print(f"C14 rate is {c14rate} Hz with Ethr = {Ethr:.2f}MeV.")


if out:
    binningE = array('d', [])
    Ewidth1, Ewidth2 = 0.01, 0.05
    for ip in range(100):
        binningE.append( Ewidth1 * ip)
    for ip in range(80):
        binningE.append( 1 + Ewidth2 * ip)
    binningE.append(5.0)
    print(binningE)
    hout = ROOT.TH2D("h1", "pES combined c14 PDF", Ntime, -0.03, 0.07, Nenergy, binningE)
    outfilename = "./Garching82703_nuePDF_NO_10kpc_pES_nuMass0.0eV_TEobs2dPDF_rebin_c14low.root"
    outfile = ROOT.TFile(outfilename, "recreate")
    
    NsiginROI, NbkginROI = 0, 0

    for it in range(pESHist.GetNbinsX()):
        t = pESHist.GetXaxis().GetBinCenter(it+1)
        twidth = pESHist.GetXaxis().GetBinWidth(it+1)
        Nsig, Nbkg = 0, 0
        for iE in range(0, pESHist.GetNbinsY(), 1):
            E = pESHist.GetYaxis().GetBinCenter(iE+1)
            Ewidth = pESHist.GetYaxis().GetBinWidth(iE+1)
            Nsig = pESHist.GetBinContent(it+1, iE+1)  # N per sec per MeV
            ## calculating C14 event rate:
            if E <= 0.2:
                Elow, Ehig = E - 0.5 * Ewidth, E + 0.5 * Ewidth
                binlow, binhig = c14Hist.FindBin(Elow), c14Hist.FindBin(Ehig)
                ratio = c14Hist.Integral(binlow, binhig) / itot
                Nbkg = ratio * nC14persec / Ewidth   # N per sec per MeV

            else:
                Nbkg = 0

            if t >= -0.02 and t <= 0.02 and E >= Ethr:
                NsiginROI += Nsig * Ewidth * twidth
                NbkginROI += Nbkg * Ewidth * twidth

            hout.SetBinContent(it+1, iE+1, Nsig+Nbkg)

    print(f"Total {NsiginROI} signals in ROI [-20, 20]ms with {Ethr:.2f}MeV.")
    print(f"Total {NbkginROI} backgrounds in ROI [-20, 20]ms with {Ethr:.2f}MeV.")
    print(f"X-bin # {hout.GetNbinsX()}, Y-bin # {hout.GetNbinsY()}.")

    hout.Write()
    outfile.Close()


















