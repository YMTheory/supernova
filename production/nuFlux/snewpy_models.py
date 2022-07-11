from astropy import units as u
from neutrino import Flavor

from astropy import units as u
from glob import glob

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import ROOT

from Fornax_2021_class import Fornax_2021

import sys

if __name__ == "__main__" :
    
    mass = 12
    nu = "nue"
    for i in range(len(sys.argv)-1):
        if sys.argv[i] == "-mass" :
            mass = int(sys.argv[i+1])
        if sys.argv[i] == "-nu":
            nu = sys.argv[i+1]

    mpl.rc('font', size=16)

    model = Fornax_2021('../../models/SNEWPY_models/Fornax_2021/lum_spec_%dM_r10000_dat.h5'%(mass))

    times = np.arange(0.00, 0.08, 0.001) * u.s
    E = np.arange(0, 101, 1) * u.MeV

    fn = "../../models/Fornax_2021/%s_%dM.dat"%(nu, mass)
    with open(fn, "w") as f:
        for t in times:
            spectra = model.get_initial_spectra(t, E)
            if nu == "nue":
                spec = spectra[Flavor.NU_E]
            elif nu == "nuebar":
                spec = spectra[Flavor.NU_E_BAR]
            elif nu == "nux":
                spec = spectra[Flavor.NU_X]
            else:
                print("Error! Wrong neutrino type !!!")
            f.write(str(t.value) + ' ')
            for j in spec:
                f.write(str(j.value) + ' ')
            f.write("\n")


    #times = np.arange(0.0, 0.06, 0.002) * u.s
    #E = np.arange(0, 101, 1) * u.MeV
    #
    #fig, axes = plt.subplots(6,5, figsize=(15,12), sharex=True, sharey=True, tight_layout=True)
    #
    #linestyles = ['-', '--', '-.', ':']
    #
    #for t, ax in zip(times, axes.flatten()):
    #    spectra = model.get_initial_spectra(t, E)
    #    for line, flavor in zip(linestyles, Flavor):
    #        if flavor == Flavor.NU_E:
    #            continue
    #        ax.plot(E, spectra[flavor], lw=3, ls=line, label=flavor.to_tex())
    #    ax.set(xlim=(0,100))
    #    ax.set_title('$t$ = {:g}'.format(t), fontsize=16)
    #    #ax.legend(loc='upper right', ncol=2, fontsize=12)
    #    ax.semilogy()
    #
    #fig.text(0.5, 0., 'energy [MeV]', ha='center')
    #fig.text(0., 0.5, f'flux [{spectra[Flavor.NU_E].unit}]', va='center', rotation='vertical');
    #plt.show()

    #### 2D neutrino E-T distribution
    #h2_nue = ROOT.TH2D("h2_nue", "visible energy v.s. T", 60, 0.0, 0.06, 100, 0, 100)
    #h2_nua = ROOT.TH2D("h2_nua", "visible energy v.s. T", 60, 0.0, 0.06, 100, 0, 100)
    #h2_nux = ROOT.TH2D("h2_nux", "visible energy v.s. T", 60, 0.0, 0.06, 100, 0, 100)

    #times = np.arange(0.00, 0.06, 0.001) * u.s
    #E = np.arange(0, 101, 1) * u.MeV
    #
    #
    #for t in times:
    #    spectra = model.get_initial_spectra(t, E)
    #    for line, flavor in zip(linestyles, Flavor):
    #        if flavor == Flavor.NU_E:
    #            for k in range(len(E)):
    #                h2_nue.SetBinContent(int((t.value)/0.001), int(E[k].value), spectra[flavor][k].value/1e50)
    #                print(t.value, int((t.value)/0.001), E[k].value, int(E[k].value), spectra[flavor][k].value/1e50)
    #        if flavor == Flavor.NU_X:
    #            for k in range(len(E)):
    #                h2_nux.SetBinContent(int((t.value)/0.001), int(E[k].value), spectra[flavor][k].value/1e50)
    #        if flavor == Flavor.NU_E_BAR:
    #            for k in range(len(E)):
    #                h2_nua.SetBinContent(int((t.value)/0.001), int(E[k].value), spectra[flavor][k].value/1e50)

    #fout = ROOT.TFile("EnuT_Burrows2D.root", "recreate")
    #h2_nue.Write()
    #h2_nua.Write()
    #h2_nux.Write()
    #fout.Close()
