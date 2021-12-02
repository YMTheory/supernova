import numpy as np
import matplotlib.pyplot as plt
import ROOT
from array import array
from math import exp
import uproot as up

def customizedPdf(tarr, parr):

    t = tarr[0]

    A    = parr[0]
    c    = parr[1]
    dT   = parr[2]
    taur = parr[3]
    taua = parr[4]

    fr = exp(-(t-dT)/taur)
    phi0 = A * exp(-(t-dT)/taua) + c

    return (1-fr) * phi0

def component1(t, dT, taur):

    fr = exp(-(t-dT)/taur)

    return (1-fr) 


def component2(t, dT, A, c, taua):
    
    phi0 = A * exp(-(t-dT)/taua) + c

    return phi0



if __name__ == "__main__" :

    #model = ROOT.TF1("model", customizedPdf, 0., 0.1, 5)
    model = ROOT.TF1("model", "(1-exp(-(x-[1])/[2])) * ([0]*exp(-(x-[1])/[3]) + [4])", 0, 0.1)

    hTime1D = ROOT.TH1D("hTime1D", "", 96, 0.004, 0.1)

    ###################################################
    ################### NO Dataset 

    model.SetParameter(0, 4.632)            # A
    model.SetParameter(1, 3.5e-3)           # dT
    model.SetParameter(2, 0.0116)           # taur
    model.SetParameter(3, 0.002)            # taub
    model.SetParameter(4, 6.71)             # c

    filename = "/junofs/users/miaoyu/supernova/wenlj/dataset/fineSpec/2.0kpc/TEvisDATA_mod82503_cha1nue_nuType-1_mh1_mNu0.0eV_2.0kpc_0.1s_Evmax25_Ethr0.1MeV_group1.root"
    inputTree = up.open(filename)["tFixedStat"]

    m_evtID = inputTree["evtID"].array()
    m_nuTime1D = inputTree["nuTime1D"].array()

    taur_arr_mh1, taub_arr_mh1 = [], []
    
    for ii in range(1, 101, 1):
        #tmp_nuTime1D = []
        print("**********************************************************")
        print("*******************Read %d SN Event*******************" %ii)
        hTime1D.Reset()
        for i, j in zip(m_evtID, m_nuTime1D):
            if i == ii:
                hTime1D.Fill(j)
                #tmp_nuTime1D.append(j)
            if i > ii:
                break
        hTime1D.Fit(model, "RE")

        taur_arr_mh1.append(model.GetParameter(3))
        taub_arr_mh1.append(model.GetParameter(4))

        print("**********************************************************")
        print(" ")
    
    ###################################################
    ################### IO Dataset 

    model.SetParameter(0, 3.136)            # A
    model.SetParameter(1, 3.8e-3)           # dT
    model.SetParameter(2, 0.0030)           # taur
    model.SetParameter(3, 0.002)            # taub
    model.SetParameter(4, 6.99)             # c


    filename = "/junofs/users/miaoyu/supernova/wenlj/dataset/fineSpec/2.0kpc/TEvisDATA_mod82503_cha1nue_nuType-1_mh2_mNu0.0eV_2.0kpc_0.1s_Evmax25_Ethr0.1MeV_group1.root"
    inputTree = up.open(filename)["tFixedStat"]

    m_evtID = inputTree["evtID"].array()
    m_nuTime1D = inputTree["nuTime1D"].array()

    taur_arr_mh2, taub_arr_mh2 = [], []
    
    for ii in range(1, 101, 1):
        #tmp_nuTime1D = []
        print("**********************************************************")
        print("*******************Read %d SN Event*******************" %ii)
        hTime1D.Reset()
        for i, j in zip(m_evtID, m_nuTime1D):
            if i == ii:
                hTime1D.Fill(j)
                #tmp_nuTime1D.append(j)
            if i > ii:
                break
        hTime1D.Fit(model, "RE")

        taur_arr_mh2.append(model.GetParameter(3))
        taub_arr_mh2.append(model.GetParameter(4))

        print("**********************************************************")
        print(" ")
    



    #for i, j in zip(m_evtID, m_nuTime1D):
    #    if i <= 1:
    #        hTime1D.Fill(j)
    #        tmp_nuTime1D.append(j)

    #hTime1D.Fit(model)

    fig, ax = plt.subplots()

    ax.hist(taur_arr_mh1, bins=20, histtype="step", lw=2, color="blue", label=r"NO $\tau_r$")
    ax.hist(taur_arr_mh2, bins=20, histtype="step", lw=2, color="blue", label=r"IO $\tau_r$")
    ax.set_xlabel("time [s]", fontsize=15)
    ax.set_ylabel("A.U.", fontsize=15)
    ax.legend(prop={"size":15})
    ax.tick_params(axis="both", which="major", labelsize=14, labelcolor="black")
    #ax.set_ylim(0, 1500)

    ax.semilogx()
    ax.semilogy()

    plt.tight_layout()
    plt.show()




    """
    func, comp1, comp2 = [], [], []
    t = np.arange(0.004, 0.101, 0.001)
    for i in t:
        func.append(model.Eval(i))
        comp1.append(component1(i, model.GetParameter(1), model.GetParameter(2)))
        comp2.append(component2(i, model.GetParameter(1), model.GetParameter(0), model.GetParameter(4), model.GetParameter(3)))



    fig, ax = plt.subplots()

    ax.hist(tmp_nuTime1D, bins=96, range=(0.004, 0.1), histtype="step", lw=2, color="blue", label="SN Data")
    ax.plot(t, func, "--", color="crimson", lw=2, label="Fitting")
    #ax.plot(t, comp1, "--", color="black", lw=2, label="Component 1")
    #ax.plot(t, comp2, "-.", color="gray", lw=2, label="Component 2")

    ax.set_xlabel("time [s]", fontsize=15)
    ax.set_ylabel("A.U.", fontsize=15)
    ax.legend(prop={"size":15})
    ax.tick_params(axis="both", which="major", labelsize=14, labelcolor="black")
    #ax.set_ylim(0, 1500)

    #ax.semilogy()

    plt.tight_layout()
    plt.savefig("primaryFit_2kpc100_mh1.pdf")
    plt.show()


    """
    








