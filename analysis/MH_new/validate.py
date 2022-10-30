#################################################
# This script save some basic distributions 
# from the fitter into pdf files for 
# validation.
# Author: Miao Yu
#################################################


import numpy as np
import ROOT
import pandas as pd


## baseline sensitivity:
def read(model, cha, E, T):

    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO = [], [], [], []
    sens_dataNO, sens_dataIO = [], []
    for i in range(10):
        df1 = pd.read_csv(f"/junofs/users/miaoyu/supernova/analysis/MH_new/results/{model}_10kpc_NO_{cha}_{E:.2f}MeV_fitTmax{T}ms_fileNo{i}.csv")
        df2 = pd.read_csv(f"/junofs/users/miaoyu/supernova/analysis/MH_new/results/{model}_10kpc_IO_{cha}_{E:.2f}MeV_fitTmax{T}ms_fileNo{i}.csv")

        tmp_deltaT_dataNO_pdfNO = df1["TbestNO"]
        tmp_deltaT_dataNO_pdfIO = df1["TbestIO"]
        tmp_sens_dataNO = df1["sens"]

        tmp_deltaT_dataIO_pdfNO = df2["TbestNO"]
        tmp_deltaT_dataIO_pdfIO = df2["TbestIO"]
        tmp_sens_dataIO = df2["sens"]

        for a, b, c, d, e, f in zip(tmp_deltaT_dataNO_pdfNO, tmp_deltaT_dataNO_pdfIO, tmp_sens_dataNO, tmp_deltaT_dataIO_pdfNO, tmp_deltaT_dataIO_pdfIO, tmp_sens_dataIO):
            deltaT_dataNO_pdfNO.append(a)
            deltaT_dataNO_pdfIO.append(b)
            sens_dataNO.append(c)
            deltaT_dataIO_pdfNO.append(d)
            deltaT_dataIO_pdfIO.append(e)
            sens_dataIO.append(f)

    deltaT_dataNO_pdfNO = np.array(deltaT_dataNO_pdfNO)
    deltaT_dataNO_pdfIO = np.array(deltaT_dataNO_pdfIO)
    sens_dataNO = np.array(sens_dataNO)
    deltaT_dataIO_pdfNO = np.array(deltaT_dataIO_pdfNO)
    deltaT_dataIO_pdfIO = np.array(deltaT_dataIO_pdfIO)
    sens_dataIO = np.array(sens_dataIO)

    return deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO 


if __name__ == "__main__":

    name = "Garching"
    model = [82503]
    chas = ["IBD", "eESIBD", "pESIBD", "pESeESIBD"]
    Tmax = np.arange(40, 61, 1)

    outfile = "./results/fitRes_validation_Garching82503.pdf"
    can = ROOT.TCanvas("can", "can", 1200, 800)
    can.Print(outfile + "[")

    for imod in model:
        imod = "Garching" + str(imod)
        for cha in chas:
            for E in [0.15]:
                for T in Tmax:
                    print(imod, cha, E, T)
                    deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO = read(imod, cha, E, T)

                    bin1, bin2 = np.min(deltaT_dataNO_pdfNO), np.max(deltaT_dataNO_pdfNO)
                    hdeltaT_dataNO_pdfNO = ROOT.TH1D("hdeltaT_dataNO_pdfNO", "", 100, bin1-1, bin2+1)
                    hdeltaT_dataNO_pdfNO.SetTitle(f"{imod}_{cha}_{E}_{T}_NOdataNOPdf")

                    bin1, bin2 = np.min(deltaT_dataNO_pdfIO), np.max(deltaT_dataNO_pdfIO)
                    hdeltaT_dataNO_pdfIO = ROOT.TH1D("hdeltaT_dataNO_pdfIO", "", 100, bin1-1, bin2+1)
                    hdeltaT_dataNO_pdfIO.SetTitle(f"{imod}_{cha}_{E}_{T}_NOdataIOPdf")
                    
                    bin1, bin2 = np.min(sens_dataNO), np.max(sens_dataNO)
                    hsens_dataNO         = ROOT.TH1D("hsens_dataNO", "", 100, bin1-1, bin2+1)
                    hsens_dataNO.SetTitle(f"{imod}_{cha}_{E}_{T}_NOdata")
                    
                    bin1, bin2 = np.min(deltaT_dataIO_pdfNO), np.max(deltaT_dataIO_pdfNO)
                    hdeltaT_dataIO_pdfNO = ROOT.TH1D("hdeltaT_dataIO_pdfNO", "", 100, bin1-1, bin2+1)
                    hdeltaT_dataIO_pdfNO.SetTitle(f"{imod}_{cha}_{E}_{T}_IOdataNOPdf")
                    
                    bin1, bin2 = np.min(deltaT_dataIO_pdfIO), np.max(deltaT_dataIO_pdfIO)
                    hdeltaT_dataIO_pdfIO = ROOT.TH1D("hdeltaT_dataIO_pdfIO", "", 100, bin1-1, bin2+1)
                    hdeltaT_dataNO_pdfIO.SetTitle(f"{imod}_{cha}_{E}_{T}_IOdataIOPdf")
                    
                    bin1, bin2 = np.min(sens_dataIO), np.max(sens_dataIO)
                    hsens_dataIO         = ROOT.TH1D("hsens_dataIO", "", 100, bin1-1, bin2+1)
                    hsens_dataIO.SetTitle(f"{imod}_{cha}_{E}_{T}_IOdata")
                    
                    for t1, t2, c1, t3, t4, c2 in zip(deltaT_dataNO_pdfNO, deltaT_dataNO_pdfIO, sens_dataNO, deltaT_dataIO_pdfNO, deltaT_dataIO_pdfIO, sens_dataIO):
                        hdeltaT_dataNO_pdfNO.Fill(t1)
                        hdeltaT_dataNO_pdfIO.Fill(t2)
                        hsens_dataNO.Fill(c1)
                        hdeltaT_dataIO_pdfNO.Fill(t3)
                        hdeltaT_dataIO_pdfIO.Fill(t4)
                        hsens_dataIO.Fill(c2)
                    
                    can.cd()
                    hdeltaT_dataNO_pdfNO.Draw()
                    can.Print(outfile)
                    can.cd()
                    hdeltaT_dataNO_pdfIO.Draw()
                    can.Print(outfile)
                    can.cd()
                    hsens_dataNO.Draw()
                    can.Print(outfile)
                    can.cd()
                    hdeltaT_dataIO_pdfNO.Draw()
                    can.Print(outfile)
                    can.cd()
                    hdeltaT_dataIO_pdfIO.Draw()
                    can.Print(outfile)
                    can.cd()
                    hsens_dataIO.Draw()
                    can.Print(outfile)
                    
    can.Print(outfile + "]")
                    
                    
                    
                    
