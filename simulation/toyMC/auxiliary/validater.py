import numpy as np
import matplotlib.pyplot as plt
import ROOT
import uproot as up

def integral2D_energy(filename, cha="eES", Ethr=0.10):
    f = ROOT.TFile(filename, "read")
    h = f.Get("h1")
    arr = np.zeros(h.GetNbinsX())
    idx0 = 0
    if cha == "pES":
        idx0 = int(Ethr / h.GetYaxis().GetBinWidth(1))
    else:
        idx0 = int( 0.20/ h.GetYaxis().GetBinWidth(1))
    for i in range(h.GetNbinsX()):
        for j in range(idx0, h.GetNbinsY(), 1):
            arr[i] += h.GetBinContent(i+1, j+1) * h.GetYaxis().GetBinWidth(1)
            
    return arr

def comparePDF(cha, MO):
    pdf1D_filename = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/Garching82703_PDF_{MO}_10kpc_{cha}_0.10MeV_newshortPDF_v2.root"
    f1d = ROOT.TFile(pdf1D_filename, "read")
    h1d = f1d.Get("h1")
    x1d, y1d = [], []
    for i in range(h1d.GetNbinsX()):
        x1d.append(h1d.GetBinCenter(i+1)*1000.)
        y1d.append(h1d.GetBinContent(i+1)/1000.)

    pdf2D_filename = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/Garching82703_PDF_{cha}_{MO}_10kpc_nuMass0.0_scale1.000_2Droot.root"
    f2d = ROOT.TFile(pdf2D_filename, "read")
    h2d = f2d.Get("h1")
    h2d1 = h2d.ProjectionX()
    x2d, y2d = [], []
    for i in range(h2d1.GetNbinsX()):
        x2d.append(h2d1.GetBinCenter(i+1)*1000)
        y2d.append(h2d1.GetBinContent(i+1)/1000.)
    y2d1 = integral2D_energy(pdf2D_filename, cha=cha, Ethr=0.10) / 1000.

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x1d, y1d, "-", label="1D")
    ax.plot(x2d, y2d1, ":", label="2D")
    ax.set_xlabel("Post-bounce time [ms]", fontsize=16)
    ax.set_ylabel("counts per ms", fontsize=16)
    ax.tick_params(axis="both", labelsize=15)
    plt.tight_layout()
    plt.show()


def compare_stat(cha):
    pdf1D_filename = f"/junofs/users/miaoyu/supernova/simulation/toyMC/scale1_poisson/Garching82703_{cha}_binneddata_NO_10kpc_thr0.10MeV_Tmin-20msTmax20ms_binning_new.root"
    pdf2D_filename = f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/Garching82703_pES_unbinneddata_NO_10.0kpc_thr0.10MeV_Tmin-20msTmax20ms_2D.root"

    f1d = up.open(pdf1D_filename)
    arr1d = f1d["binned"]["TbinConts"].array()
    NperEvt1d = []
    for subarr in arr1d:
        NperEvt1d.append(len(subarr))
    Nevt1d_mean = np.mean(NperEvt1d)

    f2d = up.open(pdf2D_filename)
    arr2d = f2d["binned"]["TbinConts"].array()
    NperEvt2d = []
    for subarr in arr2d:
        NperEvt2d.append(len(subarr))
    Nevt2d_mean = np.mean(NperEvt2d)

    print(f"====> One dimensional PDF statistics: {Nevt1d_mean}")
    print(f"====> Two dimensional PDF statistics: {Nevt2d_mean}")



if __name__ == "__main__":
    comparePDF("pES", "IO")
    #compare_stat("IBD")











