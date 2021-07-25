import numpy as np
import matplotlib.pyplot as plt
import SNGarchingIntegFcn as gar
from SNnumGarchingSrc import SNnumGarchingSrc
import SNnueXS as nuexs

imode = 82503
gar.readFluxGraph(imode)
nuTypeName = [r"$\nu_e$", r"$\nu_{\bar{e}}$", r"$\nu_x$", r"$\nu_{\bar{x}}$"]

def Lum_vs_time_allType():
    time_nue_burst, lum_nue_burst = [], []
    time_nue_acc, lum_nue_acc = [], []
    time_nue_cool, lum_nue_cool = [], []
    time_nuebar_burst, lum_nuebar_burst = [], []
    time_nuebar_acc, lum_nuebar_acc = [], []
    time_nuebar_cool, lum_nuebar_cool = [], []
    time_nux_burst, lum_nux_burst = [], []
    time_nux_acc, lum_nux_acc = [], []
    time_nux_cool, lum_nux_cool = [], []

    for i in range(gar.grLuminosity[0].GetN()):
        t = gar.grLuminosity[0].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nue_burst.append(gar.grLuminosity[0].GetPointX(i))
            lum_nue_burst.append(gar.grLuminosity[0].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nue_acc.append(gar.grLuminosity[0].GetPointX(i))
            lum_nue_acc.append(gar.grLuminosity[0].GetPointY(i))
        if 0.8 < t < 8:
            time_nue_cool.append(gar.grLuminosity[0].GetPointX(i))
            lum_nue_cool.append(gar.grLuminosity[0].GetPointY(i))

    for i in range(gar.grLuminosity[1].GetN()):
        t = gar.grLuminosity[1].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nuebar_burst.append(gar.grLuminosity[1].GetPointX(i))
            lum_nuebar_burst.append(gar.grLuminosity[1].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nuebar_acc.append(gar.grLuminosity[1].GetPointX(i))
            lum_nuebar_acc.append(gar.grLuminosity[1].GetPointY(i))
        if 0.8 < t < 8:
            time_nuebar_cool.append(gar.grLuminosity[1].GetPointX(i))
            lum_nuebar_cool.append(gar.grLuminosity[1].GetPointY(i))


    for i in range(gar.grLuminosity[2].GetN()):
        t = gar.grLuminosity[2].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nux_burst.append(gar.grLuminosity[2].GetPointX(i))
            lum_nux_burst.append(gar.grLuminosity[2].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nux_acc.append(gar.grLuminosity[2].GetPointX(i))
            lum_nux_acc.append(gar.grLuminosity[2].GetPointY(i))
        if 0.8 < t < 8:
            time_nux_cool.append(gar.grLuminosity[2].GetPointX(i))
            lum_nux_cool.append(gar.grLuminosity[2].GetPointY(i))

    plt.figure(0)
    plt.plot(time_nue_burst, lum_nue_burst, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_burst, lum_nuebar_burst, "-", label=nuTypeName[1])
    plt.plot(time_nux_burst, lum_nux_burst, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"Luminosity $10^{51}$ erg/s")
    plt.savefig("./outputs/Lum_vs_time_allType_burst_Garching%d.pdf"%imode)

    plt.figure(1)
    plt.plot(time_nue_acc, lum_nue_acc, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_acc, lum_nuebar_acc, "-", label=nuTypeName[1])
    plt.plot(time_nux_acc, lum_nux_acc, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"Luminosity $10^{51}$ erg/s")
    plt.savefig("./outputs/Lum_vs_time_allType_accretion_Garching%d.pdf"%imode)

    plt.figure(2)
    plt.plot(time_nue_cool, lum_nue_cool, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_cool, lum_nuebar_cool, "-", label=nuTypeName[1])
    plt.plot(time_nux_cool, lum_nux_cool, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel(r"Luminosity $10^{51}$ erg/s")
    plt.savefig("./outputs/Lum_vs_time_allType_cooling_Garching%d.pdf" %imode)

    print("Lum_vs_time_allType.pdf has been created !!!")


def AveE_vs_time_allType():

    time_nue_burst, aveE_nue_burst = [], []
    time_nue_acc, aveE_nue_acc = [], []
    time_nue_cool, aveE_nue_cool = [], []
    time_nuebar_burst, aveE_nuebar_burst = [], []
    time_nuebar_acc, aveE_nuebar_acc = [], []
    time_nuebar_cool, aveE_nuebar_cool = [], []
    time_nux_burst, aveE_nux_burst = [], []
    time_nux_acc, aveE_nux_acc = [], []
    time_nux_cool, aveE_nux_cool = [], []

    for i in range(gar.grAverageE[0].GetN()):
        t = gar.grAverageE[0].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nue_burst.append(gar.grAverageE[0].GetPointX(i))
            aveE_nue_burst.append(gar.grAverageE[0].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nue_acc.append(gar.grAverageE[0].GetPointX(i))
            aveE_nue_acc.append(gar.grAverageE[0].GetPointY(i))
        if 0.8 < t < 8:
            time_nue_cool.append(gar.grAverageE[0].GetPointX(i))
            aveE_nue_cool.append(gar.grAverageE[0].GetPointY(i))

    for i in range(gar.grAverageE[1].GetN()):
        t = gar.grAverageE[1].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nuebar_burst.append(gar.grAverageE[1].GetPointX(i))
            aveE_nuebar_burst.append(gar.grAverageE[1].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nuebar_acc.append(gar.grAverageE[1].GetPointX(i))
            aveE_nuebar_acc.append(gar.grAverageE[1].GetPointY(i))
        if 0.8 < t < 8:
            time_nuebar_cool.append(gar.grAverageE[1].GetPointX(i))
            aveE_nuebar_cool.append(gar.grAverageE[1].GetPointY(i))


    for i in range(gar.grAverageE[2].GetN()):
        t = gar.grAverageE[2].GetPointX(i)
        if -0.02 < t < 0.06:
            time_nux_burst.append(gar.grAverageE[2].GetPointX(i))
            aveE_nux_burst.append(gar.grAverageE[2].GetPointY(i))
        if 0.06 < t < 0.8:
            time_nux_acc.append(gar.grAverageE[2].GetPointX(i))
            aveE_nux_acc.append(gar.grAverageE[2].GetPointY(i))
        if 0.8 < t < 8:
            time_nux_cool.append(gar.grAverageE[2].GetPointX(i))
            aveE_nux_cool.append(gar.grAverageE[2].GetPointY(i))

    plt.figure(0)
    plt.plot(time_nue_burst, aveE_nue_burst, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_burst, aveE_nuebar_burst, "-", label=nuTypeName[1])
    plt.plot(time_nux_burst, aveE_nux_burst, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("AverageE/MeV")
    plt.ylim(6, 18)
    plt.savefig("./outputs/AveE_vs_time_allType_burst_Garching%d.pdf" %imode)

    plt.figure(1)
    plt.plot(time_nue_acc, aveE_nue_acc, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_acc, aveE_nuebar_acc, "-", label=nuTypeName[1])
    plt.plot(time_nux_acc, aveE_nux_acc, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("AverageE/MeV")
    plt.ylim(6, 18)
    plt.savefig("./outputs/AveE_vs_time_allType_accretion_Garching%d.pdf"%imode)

    plt.figure(2)
    plt.plot(time_nue_cool, aveE_nue_cool, "-", label=nuTypeName[0])
    plt.plot(time_nuebar_cool, aveE_nuebar_cool, "-", label=nuTypeName[1])
    plt.plot(time_nux_cool, aveE_nux_cool, "-", label=nuTypeName[2])
    plt.legend()
    plt.xlabel("time/s")
    plt.ylabel("AverageE/MeV")
    plt.ylim(6, 18)
    plt.savefig("./outputs/AveE_vs_time_allType_cooling_Garching%d.pdf"%imode)

    print("AveE_vs_time_allType.pdf has been created !!!")




def Energy_spectra_allType():
    SNGar = SNnumGarchingSrc(imode, 10)
    SNGar.readSpectrum()
    E_nue, fluence_nue, fluence_nue_NO, fluence_nue_IO              = [], [], [], []
    E_nuebar, fluence_nuebar, fluence_nuebar_NO, fluence_nuebar_IO  = [], [], [], []
    E_nux, fluence_nux, fluence_nux_NO, fluence_nux_IO              = [], [], [], []

    E_nue = np.arange(0.1, 40, 0.1)
    for i in E_nue:
        fluence_nue.append(SNGar.oneSNFluenceDet(i, 0))
        fluence_nue_NO.append(SNGar.oneSNFluenceDetMH(i, 0, 1))
        fluence_nue_IO.append(SNGar.oneSNFluenceDetMH(i, 0, 2))

    E_nuebar = np.arange(0.1, 40, 0.1)
    for i in E_nuebar:
        fluence_nuebar.append(SNGar.oneSNFluenceDet(i, 1))
        fluence_nuebar_NO.append(SNGar.oneSNFluenceDetMH(i, 1, 1))
        fluence_nuebar_IO.append(SNGar.oneSNFluenceDetMH(i, 1, 2))


    E_nux = np.arange(0.1, 40, 0.1)
    for i in E_nux:
        fluence_nux.append(SNGar.oneSNFluenceDet(i, 2))
        fluence_nux_NO.append(SNGar.oneSNFluenceDetMH(i, 2, 1))
        fluence_nux_IO.append(SNGar.oneSNFluenceDetMH(i, 2, 2))


    plt.plot(E_nue, fluence_nue,    "-" , color="blue" , label=nuTypeName[0])
    plt.plot(E_nue, fluence_nue_NO, "--", color="blue")
    plt.plot(E_nue, fluence_nue_IO, "-.", color="blue")
    plt.plot(E_nuebar, fluence_nuebar,    "-" , color="seagreen", label=nuTypeName[1])
    plt.plot(E_nuebar, fluence_nuebar_NO, "--", color="seagreen")
    plt.plot(E_nuebar, fluence_nuebar_IO, "-.", color="seagreen")
    plt.plot(E_nux, fluence_nux,    "-" ,color="orange", label=nuTypeName[2])
    plt.plot(E_nux, fluence_nux_NO, "--",color="orange")
    plt.plot(E_nux, fluence_nux_IO, "-.",color="orange")

    plt.hlines(1.35e10, 20, 22, linestyle="-")
    plt.text(23, 1.35e10, "w/o osc")
    plt.hlines(1.25e10, 20, 22, linestyle="--")
    plt.text(23, 1.25e10, "w/ osc, NO")
    plt.hlines(1.15e10, 20, 22, linestyle="-.")
    plt.text(23, 1.15e10, "w/ osc, IO")

    plt.xlabel(r"$E_\nu$/MeV")
    plt.ylabel("fluence")
    plt.legend()
    plt.savefig("./outputs/Energy_spectra_allType_Garching%d.pdf"%imode)




import ROOT
from array import array
def Energy_time_2D_OneType(tp):
    # time binning
    tmin, tmax = -0.1, 0.1
    #t1_thr, t2_thr = -0.025, 0.05
    t1_thr, t2_thr = -0.025, 0.1
    step_t0, step_t1, step_t2 = 0.00005, 0.00005, 0.0005
    npt0 = round( (t1_thr-tmin)/step_t0 )
    npt1 = round( (t2_thr-t1_thr)/step_t1 )
    npt2 = round( (tmax-t2_thr)/step_t2 )
    nbin_time = npt0 + npt1 + npt2
    
    # time binning
    binning_time = array('d', [])
    for ip in range(npt0):
        binning_time.append( tmin+step_t0*ip )
    
    for ip in range(npt0, npt0 + npt1):
        binning_time.append( t1_thr+step_t1*(ip-npt0) )

    for ip in range(npt0 + npt1, npt0 + npt1 + npt2):
        binning_time.append( t2_thr+step_t2*(ip-npt0-npt1) )
    
    binning_time.append(tmax)

    # Enu binning
    Evmin,   Evmax,   step_Ev   = 0.0, 60.0, 0.2
    nbins_Ev   = round( (Evmax-Evmin)/step_Ev )


    hist = ROOT.TH2D("hET_source", "", nbin_time, binning_time, nbins_Ev, Evmin, Evmax)
    for i in range(hist.GetNbinsX()):
        for j in range(hist.GetNbinsY()):
            t = hist.GetXaxis().GetBinCenter(i+1)
            E = hist.GetYaxis().GetBinCenter(j+1)
            cont = gar.getEventAtTime(t, E, tp)
            hist.Fill(t, E, cont)

    fout = ROOT.TFile("TEv_source_Garching%d.root"%imode, "recreate")
    hist.Write()
    fout.Close()
    print("TEv_source_Garching rootfile has been created !!!")


def totalXS():
    totxs0, totxs1, totxs2, totxs3 = [], [], [], []
    Ev = np.arange(0.0, 60, 0.1)
    for i in Ev:
        totxs0.append( nuexs.totalXS(i, 0) )
        totxs1.append( nuexs.totalXS(i, 1) )
        totxs2.append( nuexs.totalXS(i, 2) )
        totxs3.append( nuexs.totalXS(i, 3) )

    plt.plot(Ev, totxs0, label=nuTypeName[0])
    plt.plot(Ev, totxs1, label=nuTypeName[1])
    plt.plot(Ev, totxs2, label=nuTypeName[2])
    plt.plot(Ev, totxs3, label=nuTypeName[3])

    plt.grid(True)
    plt.xlabel(r"$T_{e^-}$/MeV")
    plt.ylabel(r"$\nu$-eES cross section[cm$^{2}$]")
    plt.legend()
    plt.savefig("./outputs/nue_totxs.pdf")



def diffXS():
    diffxs0, diffxs1, diffxs2, diffxs3 = [], [], [], []
    Enu = 15  # MeV
    Evis = np.arange(0.1, Enu-0.1, 0.1)
    for i in Evis:
        diffxs0.append( nuexs.differentialXS(Enu, i, 0) )
        diffxs1.append( nuexs.differentialXS(Enu, i, 1) )
        diffxs2.append( nuexs.differentialXS(Enu, i, 2) )
        diffxs3.append( nuexs.differentialXS(Enu, i, 3) )

    plt.plot(Evis, diffxs0, label=nuTypeName[0])
    plt.plot(Evis, diffxs1, label=nuTypeName[1])
    plt.plot(Evis, diffxs2, label=nuTypeName[2])
    plt.plot(Evis, diffxs3, label=nuTypeName[3])

    plt.grid(True)
    plt.xlabel(r"$T_{e^-}$/MeV")
    plt.ylabel(r"$\nu$-eES cross section[cm$^{2}$ MeV$^{-1}$]")
    plt.legend()
    plt.savefig("./outputs/nue15MeV_diffxs.pdf")



def NEv_spectra_allType():
    SNGar = SNnumGarchingSrc(imode, 10)
    SNGar.readSpectrum()
    T_nue, nEv_nue, nEv_nue_NO, nEv_nue_IO              = [], [], [], []
    T_nuebar, nEv_nuebar, nEv_nuebar_NO, nEv_nuebar_IO  = [], [], [], []
    T_nux, nEv_nux, nEv_nux_NO, nEv_nux_IO              = [], [], [], []

    T_nue = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nue:
        nEv_nue.append(SNGar.oneSNFluenceDetTimeIntegE(i, 0, 0))
        nEv_nue_NO.append(SNGar.oneSNFluenceDetTimeIntegE(i, 0, 1))
        nEv_nue_IO.append(SNGar.oneSNFluenceDetTimeIntegE(i, 0, 2))

    T_nuebar = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nuebar:
        nEv_nuebar.append(SNGar.oneSNFluenceDetTimeIntegE(i, 1, 0))
        nEv_nuebar_NO.append(SNGar.oneSNFluenceDetTimeIntegE(i, 1, 1))
        nEv_nuebar_IO.append(SNGar.oneSNFluenceDetTimeIntegE(i, 1, 2))


    T_nux = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nux:
        nEv_nux.append(SNGar.oneSNFluenceDetTimeIntegE(i, 2, 0))
        nEv_nux_NO.append(SNGar.oneSNFluenceDetTimeIntegE(i, 2, 1))
        nEv_nux_IO.append(SNGar.oneSNFluenceDetTimeIntegE(i, 2, 2))


    plt.plot(T_nue, nEv_nue,    "-" , color="blue" , label=nuTypeName[0]+", w/o osc")
    plt.plot(T_nue, nEv_nue_NO, "--", color="blue" , label=nuTypeName[0]+", w/ NO")
    plt.plot(T_nue, nEv_nue_IO, "-.", color="blue" , label=nuTypeName[0]+", w/ IO")
    #plt.plot(T_nuebar, nEv_nuebar,    "-" , color="seagreen", label=nuTypeName[1])
    #plt.plot(T_nuebar, nEv_nuebar_NO, "--", color="seagreen")
    #plt.plot(T_nuebar, nEv_nuebar_IO, "-.", color="seagreen")
    #plt.plot(T_nux, nEv_nux,    "-" ,color="orange", label=nuTypeName[2])
    #plt.plot(T_nux, nEv_nux_NO, "--",color="orange")
    #plt.plot(T_nux, nEv_nux_IO, "-.",color="orange")

    plt.xlabel("time/s")
    plt.ylabel("nEv")
    plt.legend()
    plt.savefig("./outputs/NEv_spectra_nue_Garching%d.pdf"%imode)



from EvEneSpec import EvEneSpec
def Energy_spectra_allType_miao():
    spec = EvEneSpec(1, 12, 1)
    Ev = np.arange(0.1, 60, 0.1)
    dFdE = spec.KRJ_dFdE(Ev)
    
    plt.plot(Ev, dFdE, "-")
    plt.show()




def compare_garching_NEvT():
    mode_arr = [71200, 73500, 81500, 82503, 92000]
    T_nue0, nEv_nue0, nEv_nue_NO0, nEv_nue_IO0              = [], [], [], []
    T_nue1, nEv_nue1, nEv_nue_NO1, nEv_nue_IO1              = [], [], [], []
    T_nue2, nEv_nue2, nEv_nue_NO2, nEv_nue_IO2              = [], [], [], []
    T_nue3, nEv_nue3, nEv_nue_NO3, nEv_nue_IO3              = [], [], [], []
    T_nue4, nEv_nue4, nEv_nue_NO4, nEv_nue_IO4              = [], [], [], []

    SNGar0 = SNnumGarchingSrc(mode_arr[0], 10)
    SNGar0.readSpectrum()
    T_nue0 = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nue0:
        nEv_nue0.append(SNGar0.oneSNFluenceDetTimeIntegE(i, 0, 0))
        nEv_nue_NO0.append(SNGar0.oneSNFluenceDetTimeIntegE(i, 0, 1))
        nEv_nue_IO0.append(SNGar0.oneSNFluenceDetTimeIntegE(i, 0, 2))

    SNGar1 = SNnumGarchingSrc(mode_arr[1], 10)
    SNGar1.readSpectrum()
    T_nue1 = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nue1:
        nEv_nue1.append(SNGar1.oneSNFluenceDetTimeIntegE(i, 0, 0))
        nEv_nue_NO1.append(SNGar1.oneSNFluenceDetTimeIntegE(i, 0, 1))
        nEv_nue_IO1.append(SNGar1.oneSNFluenceDetTimeIntegE(i, 0, 2))



    SNGar2 = SNnumGarchingSrc(mode_arr[2], 10)
    SNGar2.readSpectrum()
    T_nue2 = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nue2:
        nEv_nue2.append(SNGar2.oneSNFluenceDetTimeIntegE(i, 0, 0))
        nEv_nue_NO2.append(SNGar2.oneSNFluenceDetTimeIntegE(i, 0, 1))
        nEv_nue_IO2.append(SNGar2.oneSNFluenceDetTimeIntegE(i, 0, 2))


    SNGar3 = SNnumGarchingSrc(mode_arr[3], 10)
    SNGar3.readSpectrum()
    T_nue3 = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nue3:
        nEv_nue3.append(SNGar3.oneSNFluenceDetTimeIntegE(i, 0, 0))
        nEv_nue_NO3.append(SNGar3.oneSNFluenceDetTimeIntegE(i, 0, 1))
        nEv_nue_IO3.append(SNGar3.oneSNFluenceDetTimeIntegE(i, 0, 2))


    SNGar4 = SNnumGarchingSrc(mode_arr[4], 10)
    SNGar4.readSpectrum()
    T_nue4 = np.arange(-0.1, 0.1, 0.0001)
    for i in T_nue4:
        nEv_nue4.append(SNGar4.oneSNFluenceDetTimeIntegE(i, 0, 0))
        nEv_nue_NO4.append(SNGar4.oneSNFluenceDetTimeIntegE(i, 0, 1))
        nEv_nue_IO4.append(SNGar4.oneSNFluenceDetTimeIntegE(i, 0, 2))



    plt.plot(T_nue0, nEv_nue_NO0,    "-"  , label=r"Garching%d, $\nu_e NO$"%mode_arr[0])
    plt.plot(T_nue1, nEv_nue_NO1,    "-"  , label=r"Garching%d, $\nu_e NO$"%mode_arr[1])
    plt.plot(T_nue2, nEv_nue_NO2,    "-"  , label=r"Garching%d, $\nu_e NO$"%mode_arr[2])
    plt.plot(T_nue3, nEv_nue_NO3,    "-"  , label=r"Garching%d, $\nu_e NO$"%mode_arr[3])
    plt.plot(T_nue4, nEv_nue_NO4,    "-"  , label=r"Garching%d, $\nu_e NO$"%mode_arr[4])

    plt.xlabel("time/s")
    plt.ylabel("nEv")
    plt.legend()

    plt.savefig("./outputs/compare_Garching_nEvt2T_NO.pdf")
    plt.show()




libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
#libSNsimDir = os.getenv('SNCODEDIR') + '/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)
print('Load', libSNsimDir)
def NEvis_spectra_allType():

    snDet = ROOT.SNdetect.instance()
    modelNum = imode
    snDet.setSrcModel(modelNum)
    modelSrc = snDet.getPointerSrc()

    # -- configurations -- #
    snDet.initChannel(1)
    modelSrc.setSNDistance(10)
    Ethr = 0.2
    snDet.getPointerEffectLS().setThresholdE(Ethr)
    snDet.initFCN()

    
    NEvis, NEvis_NO, NEvis_IO = [], [], []
    T_vis = np.arange(-0.1, 0.1, 0.001)
    for i in T_vis:
        NEvis.append(snDet.getEventAboveEthrVisAtTime(i, Ethr, -1, 0))
        NEvis_NO.append(snDet.getEventAboveEthrVisAtTime(i, Ethr, -1, 1))
        NEvis_IO.append(snDet.getEventAboveEthrVisAtTime(i, Ethr, -1, 2))
        print(i, NEvis[-1], NEvis_NO[-1], NEvis_IO[-1])
    
    plt.plot(T_vis, NEvis, "-", label="w/ osc")
    plt.plot(T_vis, NEvis_NO, "--", label="NO")
    plt.plot(T_vis, NEvis_IO, "-.", label="IO")
    plt.legend()
    plt.xlabel(r"$E_{dep}$/MeV")
    plt.ylabel("Event Rate")
    plt.savefig("./outputs/NEvis_spectra_allType.pdf")
    plt.show()





