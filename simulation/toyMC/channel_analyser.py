#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Fri Mar  3 15:23:13 2023
# File Name: channel_analyser.py
"""

import uproot as up
import numpy as np
import os, sys
import uproot as up
from collections import Counter
import scipy.integrate as integrate
from scipy import interpolate
from scipy.special import gamma
#import ROOT
import time
from tqdm import tqdm

class channel :
    """
    Define a new class describing a specific detection channel in JUNO.
    """
    
    def __init__(self, name:str, MH:str, model:str, modelNo:int, Ethr:float, dist=10, fitTmin=10, fitTmax=50, fitEmax=80, fileNo=0, exp="", nuMass=0.0) -> None:
        self.name   = name
        self.MH     = MH
        self.model  = model
        self.modelNo = modelNo
        self.Ethr   = Ethr
        self.dist   = dist
        self.scale  = 10. * 10. / dist / dist
        self.nuMass = nuMass

        self.fitTmin = fitTmin
        self.fitTmax = fitTmax
        self.fitEmax = fitEmax

        self.fileNo  = fileNo
        self.NevtPerFile = 1000
        
        production_path = "/junofs/users/miaoyu/supernova/production/"
        #self.datafile   = f"{production_path}/Data/{dist}kpc/{exp}/{model}{modelNo}_{name}_data_{MH}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_new_{fileNo}.root"
        #self.datafile   = f"{production_path}/Data/{dist}kpc/{model}/{model}{modelNo}_{name}_data_{MH}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_merger.root"
        #self.datafile   = f"{production_path}/Data/{dist}kpc/{model}/{model}{modelNo}_{name}_data_{MH}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_new_{fileNo}.root"
        Ethrm = int(Ethr * 1e3)
        self.datafile   = f"/junofs/users/miaoyu/supernova/simulation/C++/Data/{name}/{Ethrm}keV/{model}{modelNo}_{name}_data_{MH}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning.root"
        # for other Exps...
        #self.pdfNOfile  = f"{production_path}/PDFs/10kpc/{exp}/{model}{modelNo}_PDF_NO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF_{exp}.root"
        ##self.pdfIOfile  = f"{production_path}/PDFs/10kpc/{exp}/{model}{modelNo}_PDF_IO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF_{exp}.root"
        #self.pdfNOfile  = f"{production_path}/PDFs/10kpc/{model}/{model}{modelNo}_PDF_NO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        #self.pdfIOfile  = f"{production_path}/PDFs/10kpc/{model}/{model}{modelNo}_PDF_IO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        #self.pdfNOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_NO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        #self.pdfIOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_IO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        #self.pdfNOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs/1D/{model}{modelNo}_PDF_NO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        #self.pdfIOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs/1D/{model}{modelNo}_PDF_IO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        
        self.pdfNOfile0 = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_NO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF_v2.root"
        self.pdfIOfile0 = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_IO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF_v2.root"

        if self.model == "Garching":
            self.pdfNOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_NO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
            self.pdfIOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_IO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        elif self.model == "Burrows":
            if self.name == "pES":
                self.pdfNOfile = f"/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Burrows/{model}{modelNo}_PDF_NO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
                self.pdfIOfile = f"/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Burrows/{model}{modelNo}_PDF_IO_10kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
            else:
                self.pdfNOfile = f"/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Burrows/{model}{modelNo}_PDF_NO_10kpc_{name}_0.20MeV_newshortPDF.root"
                self.pdfIOfile = f"/junofs/users/miaoyu/supernova/production/PDFs/10kpc/Burrows/{model}{modelNo}_PDF_IO_10kpc_{name}_0.20MeV_newshortPDF.root"


        self.data2Dfile = f"/junofs/users/miaoyu/supernova/simulation/toyMC/Data2d/{model}{modelNo}_{name}_unbinneddata_{MH}_{dist:.1f}kpc_thr{Ethr:.2f}MeV_Tmin-20msTmax{fitTmax}ms_2D.root"
        self.pdf2DNOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/{model}{modelNo}_PDF_{name}_NO_{dist}kpc_nuMass{nuMass:.1f}_scale1.000_test2Dnew.root"
        self.pdf2DIOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/{model}{modelNo}_PDF_{name}_IO_{dist}kpc_nuMass{nuMass:.1f}_scale1.000_test2Dnew.root"
        self.pdf2DNOfile0 = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/{model}{modelNo}_PDF_{name}_NO_{dist}kpc_nuMass0.0_scale1.000_test2Dnew.root"
        self.pdf2DIOfile0 = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs2d/{model}{modelNo}_PDF_{name}_IO_{dist}kpc_nuMass0.0_scale1.000_test2Dnew.root"

        ####### Datasets and PDFs
        self.data_array = None
        self.dataT_array = None
        self.dataE_array = None
        self.binned_data_array = None
        self.nentries   = 0
        self.startEvt   = 0
        self.endEvt     = self.NevtPerFile
        #self.pdfNO      = None
        #self.pdfIO      = None
        self.pdfNOx     = None
        self.pdfNOy     = None
        self.pdfIOx     = None
        self.pdfIOy     = None
        # for 1D absolute mass analysis
        self.pdfNOx0    = None
        self.pdfNOy0    = None
        self.pdfIOx0    = None
        self.pdfIOy0    = None
        # for 2D MO analysis
        self.pdf2DNOT     = None
        self.pdf2DNOE     = None
        self.pdf2DNO     = None
        self.f2dNO      = None
        self.pdf2DIOT     = None
        self.pdf2DIOE     = None
        self.pdf2DIO     = None
        self.f2dIO      = None
        # for 2D absolute mass analysis
        self.pdf2DNOT0    = None
        self.pdf2DNOE0    = None
        self.pdf2DNO0    = None
        self.f2dNO0     = None
        self.pdf2DIOT0    = None
        self.pdf2DIOE0    = None
        self.pdf2DIO0    = None
        self.f2dIO0     = None

        self.glow       = None
        self.ghig       = None
        self.c14rate    = 0

        self.Tbinwidth   = 0.00001
        self.Ebinwidth  = 0.2
        if self.name == "pES":
            self.Ebinwidth = 0.05

        self.glow = None
        self.ghig = None

    def setC14rate(self, rate):
        self.c14rate = rate

    def setNevtPerFile(self, N) -> None:
        self.NevtPerFile = N

    def getNevtPerFile(self) -> int:
        return self.NevtPerFile
    
    def setDataFilePath(self, datafile) -> None:
        self.datafile = datafile

    def setData2DFilePath(self, datafile) -> None:
        self.data2Dfile = datafile

    def setNOPdfFilePath(self, pdffile) -> None:
        self.pdfNOfile = pdffile

    def setIOPdfFilePath(self, pdffile) -> None:
        self.pdfIOfile = pdffile

    def setNOPdfFile0Path(self, pdffile) -> None:
        self.pdfNOfile0 = pdffile

    def setIOPdfFile0Path(self, pdffile) -> None:
        self.pdfIOfile0 = pdffile

    def setNOPdf2DFilePath(self, pdffile) -> None:
        self.pdf2DNOfile = pdffile

    def setIOPdf2DFilePath(self, pdffile) -> None:
        self.pdf2DIOfile = pdffile

    def setNOPdf2DFile0Path(self, pdffile) -> None:
        self.pdf2DNOfile0 = pdffile

    def setIOPdf2DFile0Path(self, pdffile) -> None:
        self.pdf2DIOfile0 = pdffile

    def setStartEvtId(self, idd):
        self.startEvt = idd
        
    def setEndEvtId(self, idd):
        self.endEvt = idd

    def _load_data_ak(self) -> None:
        """
        Load ttree ak array data from root files.
        """
        try:
            with up.open(self.datafile) as f:
                print("--------------------------------------------------------------")
                print(f"Load datafile of channel {self.name} from ->")
                print(self.datafile)
                print("--------------------------------------------------------------")
                numEntries = f["binned"]["TbinConts"].num_entries
                if self.endEvt == 0:
                    self.startEvt = 0
                    self.endEvt = numEntries # load all entries from the TBranch.
                nuTime = f["binned"]["TbinConts"].array(entry_start=self.startEvt, entry_stop=self.endEvt)
                evtNum = self.endEvt - self.startEvt
                print(f"Total loaded event number = {evtNum} from Event {self.startEvt} to {self.endEvt}.")
                self.data_array = nuTime
                self.nentries = evtNum
        except FileExistsError:
            print(f"The data file {self.datafile} dose not exist! :(")
            sys.exit(-1)


    def _load_data_binned(self) -> None:
        """
        Load binned data contents.
        """
        try:
            with up.open(self.datafile) as f:
                print("--------------------------------------------------------------")
                print(f"Load binned datafile of channel {self.name} from ->")
                print(self.datafile)
                print("--------------------------------------------------------------")
                numEntries = f["binned"]["TbinConts"].num_entries
                if self.endEvt == 0:
                    self.startEvt = 0
                    self.endEvt = numEntries # load all entries from the TBranch.
                nuTime = f["binned"]["TbinConts"].array(entry_start=self.startEvt, entry_stop=self.endEvt)
                evtNum = self.endEvt - self.startEvt
                print(f"Total loaded event number = {evtNum} from Event {self.startEvt} to {self.endEvt}.")
                self.binned_data_array = nuTime
                self.nentries = evtNum
        except FileExistsError:
            print(f"The data file {self.datafile} dose not exist! :(")
            sys.exit(-1)


    def _load_data(self) -> None:
        """
        Load ttree data from root files.
        """
        try:
            with up.open(self.datafile) as f:
                print(self.datafile)
                numEntries = f["tFixedStat"]["evtID"].num_entries
                #numEntries = f["binned"]["TbinCont"].num_entries
                nuPerEvent = int(numEntries / self.NevtPerFile)
                print(f"Statistics of channel {self.name} is {nuPerEvent}.")
                if self.endEvt == 0:
                    self.startEvt = 0
                    self.endEvt = numEntries
                #nuTime = f["binned"]["TbinCont"].array(entry_start=self.startEvt*nuPerEvent, entry_stop=self.endEvt*nuPerEvent)
                nuTime = f["tFixedStat"]["nuTime1D"].array(entry_start=self.startEvt*nuPerEvent, entry_stop=self.endEvt*nuPerEvent)
                evtNum = self.endEvt - self.startEvt
                print(f"Total loaded event number = {evtNum} from Event {self.startEvt} to {self.endEvt}.")
                nuTime = np.reshape(nuTime, (evtNum, nuPerEvent))
                nuTime = np.array(nuTime)
        except FileNotFoundError:
            print(f"The data file {self.datafile} dose not exist! :(")
            sys.exit(-1)

        self.data_array = nuTime
        self.nentries   = self.startEvt - self.endEvt

    def _load_data2D(self) -> None:
        """
        Load ttree ak array data from root files.
        """
        try:
            with up.open(self.data2Dfile) as f:
                print("--------------------------------------------------------------")
                print(f"Load datafile of channel {self.name} from ->")
                print(self.data2Dfile)
                print("--------------------------------------------------------------")
                numEntries = f["binned"]["TbinConts"].num_entries
                if self.endEvt == 0:
                    self.startEvt = 0
                    self.endEvt = numEntries # load all entries from the TBranch.
                nuTime = f["binned"]["TbinConts"].array(entry_start=self.startEvt, entry_stop=self.endEvt)
                nuEnergy = f["binned"]["EbinConts"].array(entry_start=self.startEvt, entry_stop=self.endEvt)
                evtNum = self.endEvt - self.startEvt
                Nsubarr = []
                for subarr in nuTime:
                    Nsubarr.append(len(subarr))
                Nsubarr = np.array(Nsubarr)
                Nave = np.mean(Nsubarr)
                Nsigma = np.std(Nsubarr)
                print(f"Total loaded event number = {evtNum} from Event {self.startEvt} to {self.endEvt} in channel {self.name}, where the average signal number in one event is {Nave} and sigma is {Nsigma}.")
                self.dataT_array = nuTime
                self.dataE_array = nuEnergy
                self.nentries = evtNum
        except FileExistsError:
            print(f"The data file {self.datafile} dose not exist! :(")
            sys.exit(-1)


    def _load_pdf(self) -> None:
        """
        Load TH1 PDFs from root files for both NO and IO cases.
        """
        print(f"\n ========= Load {self.name} 1D PDF ========= \n")
        try:
            print(self.pdfNOfile)
            f = up.open(self.pdfNOfile)
            tmp_h1 = f["h1"] 
        except FileNotFoundError:
            print(f"The pdf file {self.pdfNOfile} dose not exist! :(")
            sys.exit(-1)
            
        # implement a new ROOT TH1D as pdfs
        axis = tmp_h1.axis()    
        #self.pdfNO = ROOT.TH1D(f"{self.name}pdfNO", "NO time pdf", len(axis.centers()), axis.low, axis.high)
        #for j, cont in enumerate(tmp_h1.values()):
        #    self.pdfNO.SetBinContent(j+1, cont)
        self.pdfNOx = axis.centers()
        self.pdfNOy = tmp_h1.values()
        
        try:
            print(self.pdfIOfile)
            f = up.open(self.pdfIOfile)
            tmp_h1 = f["h1"] 
        except FileNotFoundError:
            print(f"The pdf file {self.pdfIOfile} dose not exist! :(")
            sys.exit(-1)
            
        # implement a new ROOT TH1D as pdfs
        axis = tmp_h1.axis()    
        #self.pdfIO = ROOT.TH1D(f"{self.name}pdfIO", "IO time pdf", len(axis.centers()), axis.low, axis.high)
        #for j, cont in enumerate(tmp_h1.values()):
        #    self.pdfIO.SetBinContent(j+1, cont)
        self.pdfIOx = axis.centers()
        self.pdfIOy = tmp_h1.values()
        

    def _load_pdf0(self) -> None:
        """
        Load TH1 PDFs from root files for both NO and IO cases.
        """
        print(f"\n ========= Load {self.name} 1D PDF0 ========= \n")
        try:
            print(self.pdfNOfile0)
            f = up.open(self.pdfNOfile0)
            tmp_h1 = f["h1"] 
        except FileNotFoundError:
            print(f"The pdf file {self.pdfNOfile0} dose not exist! :(")
            sys.exit(-1)
            
        # implement a new ROOT TH1D as pdfs
        axis = tmp_h1.axis()    
        #self.pdfNO = ROOT.TH1D(f"{self.name}pdfNO", "NO time pdf", len(axis.centers()), axis.low, axis.high)
        #for j, cont in enumerate(tmp_h1.values()):
        #    self.pdfNO.SetBinContent(j+1, cont)
        self.pdfNOx0 = axis.centers()
        self.pdfNOy0 = tmp_h1.values()
        
        try:
            print(self.pdfIOfile0)
            f = up.open(self.pdfIOfile0)
            tmp_h1 = f["h1"] 
        except FileNotFoundError:
            print(f"The pdf file {self.pdfIOfile0} dose not exist! :(")
            sys.exit(-1)
            
        # implement a new ROOT TH1D as pdfs
        axis = tmp_h1.axis()    
        #self.pdfIO = ROOT.TH1D(f"{self.name}pdfIO", "IO time pdf", len(axis.centers()), axis.low, axis.high)
        #for j, cont in enumerate(tmp_h1.values()):
        #    self.pdfIO.SetBinContent(j+1, cont)
        self.pdfIOx0 = axis.centers()
        self.pdfIOy0 = tmp_h1.values()
        


    def _load_pdf2D(self) -> None:
        print(f"\n ========= Load {self.name} 2D PDF ========= \n")
        try:
            print(self.pdf2DNOfile)
            f = up.open(self.pdf2DNOfile)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            print(f"THe pdf file {self.pdf2DNOfile} does not exist! :(")
            sys.exit(-1)
        xaxis = tmp_h1.axis("x")    
        yaxis = tmp_h1.axis("y")
        self.pdf2DNOT = xaxis.centers()
        self.pdf2DNOE = yaxis.centers()
        self.pdf2DNO  = tmp_h1.values()
        self.f2dNO = interpolate.interp2d(self.pdf2DNOT, self.pdf2DNOE, self.pdf2DNO.T, kind="linear")
             
        try:
            print(self.pdf2DIOfile)
            f = up.open(self.pdf2DIOfile)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            print(f"THe pdf file {self.pdf2DIOfile} does not exist! :(")
            sys.exit(-1)
        xaxis = tmp_h1.axis("x")    
        yaxis = tmp_h1.axis("y")
        self.pdf2DIOT = xaxis.centers()
        self.pdf2DIOE = yaxis.centers()
        self.pdf2DIO  = tmp_h1.values()
        self.f2dIO = interpolate.interp2d(self.pdf2DIOT, self.pdf2DIOE, self.pdf2DIO.T, kind="linear")
        
    
    def _load_pdf2D0(self) -> None:
        print(f"\n ========= Load {self.name} 2D PDF0 ========= \n")
        try:
            print(self.pdf2DNOfile0)
            f = up.open(self.pdf2DNOfile0)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            print(f"THe pdf file {self.pdf2DNOfile} does not exist! :(")
            sys.exit(-1)
        xaxis = tmp_h1.axis("x")    
        yaxis = tmp_h1.axis("y")
        self.pdf2DNOT0 = xaxis.centers()
        self.pdf2DNOE0 = yaxis.centers()
        self.pdf2DNO0  = tmp_h1.values()
        self.f2dNO0 = interpolate.interp2d(self.pdf2DNOT0, self.pdf2DNOE0, self.pdf2DNO0.T, kind="linear")
             
        try:
            print(self.pdf2DIOfile0)
            f = up.open(self.pdf2DIOfile0)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            print(f"THe pdf file {self.pdf2DIOfile} does not exist! :(")
            sys.exit(-1)
        xaxis = tmp_h1.axis("x")    
        yaxis = tmp_h1.axis("y")
        self.pdf2DIOT0 = xaxis.centers()
        self.pdf2DIOE0 = yaxis.centers()
        self.pdf2DIO0  = tmp_h1.values()
        self.f2dIO0 = interpolate.interp2d(self.pdf2DIOT0, self.pdf2DIOE0, self.pdf2DIO0.T, kind="linear")



        
    def get_one_event(self, event_id:int):
        """
        get a single event with ID == event_id.
        """
        event_id = event_id - self.startEvt
        return self.data_array[event_id]
    
    def get_one_binned_event(self, event_id:int):
        event_id = event_id - self.startEvt
        return self.binned_data_array[event_id]

    def get_one_event1D(self, event_id:int):
        event_id = event_id - self.startEvt
        return self.dataT_array[event_id]
    
    def get_one_event2D(self, event_id:int):
        event_id = event_id - self.startEvt
        return self.dataT_array[event_id], self.dataE_array[event_id]
    
    def _pdfNO_func0(self, t):
        return np.interp(t, self.pdfNOx0, self.pdfNOy0) 
    
    def _pdfIO_func0(self, t):
        return np.interp(t, self.pdfIOx0, self.pdfIOy0)

    
    def _pdfNO_func(self, t):
        return np.interp(t, self.pdfNOx, self.pdfNOy) 
    
    def _pdfIO_func(self, t):
        return np.interp(t, self.pdfIOx, self.pdfIOy)

    def _pdf2DNO_func(self, t, E):
        return self.f2dNO(t, E)

    def _pdf2DIO_func(self, t, E):
        return self.f2dIO(t, E) 

    def _pdf2DNO_func0(self, t, E):
        return self.f2dNO0(t, E)

    def _pdf2DIO_func0(self, t, E):
        return self.f2dIO0(t, E) 

    def getNsigCurrentEvent(self, evtid:int) -> int:
        evtid = evtid - self.startEvt
        return len(self.data_array[evtid])

    
    def calc_NLL_NO(self, data, dT) -> float:
        """
        calculate NLL for given dataset and PDF, where PDF is shifted by dT.
        """
        #time_start = time.time()
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for i in data:
            #tmp_nll = self.pdfNO.Interpolate(i + dT)
            tmp_nll = np.interp(i+dT, self.pdfNOx, self.pdfNOy) * self.scale
            if tmp_nll <= 0:
                continue
                #tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        #bin1 = self.pdfNO.GetXaxis().FindBin(tmin)
        #bin2 = self.pdfNO.GetXaxis().FindBin(tmax)
        #intg = self.pdfNO.Integral(bin1, bin2, "width")

        intg = integrate.quad(self._pdfNO_func, tmin, tmax)[0] * self.scale
        nll -= intg    
        #time_end = time.time()
        #print(f"Time consumed = {time_end - time_start} s.")
        return -nll


    def calc_NLL_IO(self, data, dT) -> float:
        """
        calculate NLL for given dataset and PDF, where PDF is shifted by dT.
        """
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for i in data:
            #tmp_nll = self.pdfIO.Interpolate(i + dT)
            tmp_nll = np.interp(i+dT, self.pdfIOx, self.pdfIOy) * self.scale
            if tmp_nll <= 0:
                #tmp_nll = 1e-10
                continue
            nll += np.log(tmp_nll)
        
        #bin1 = self.pdfIO.GetXaxis().FindBin(tmin)
        #bin2 = self.pdfIO.GetXaxis().FindBin(tmax)
        #intg = self.pdfIO.Integral(bin1, bin2, "width")
        
        intg = integrate.quad(self._pdfIO_func, tmin, tmax)[0] * self.scale
        
        nll -= intg 

        return -nll


    def calc_NLL_NO2D(self, dataT, dataE, dT) -> float: ## input data time unit: ms
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for t, E in zip(dataT, dataE):
            tmp_nll = self._pdf2DNO_func((t+dT), E)[0] * self.scale
            if tmp_nll <= 0:
                tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        intg = integrate.quad(self._pdfNO_func, tmin, tmax)[0] * self.scale
        nll -= intg
        return -nll
    
    
    def calc_NLL_IO2D(self, dataT, dataE, dT) -> float:
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for t, E in zip(dataT, dataE):
            tmp_nll = self._pdf2DIO_func((t+dT), E)[0]  * self.scale
            if tmp_nll <= 0:
                tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        intg = integrate.quad(self._pdfIO_func, tmin, tmax)[0] * self.scale
        nll -= intg
        return -nll


    def calc_binnedNLL_NO(self, data, dT) -> float:
        nll = 0
        stepT = 0.01
        # bins with counts
        t_data = np.where(data) 
        n = data[t_data]
        t_data = np.array(t_data) * stepT + self.fitTmin
        t_pdf = t_data + dT 
        s = self._pdfNO_func(t_pdf) * stepT * self.scale
        tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
        tmp_nll = np.sum(tmp_nll) 
        # bins without counts
        t_data0 = np.where(data == 0)
        t_pdf0 = self.fitTmin + np.array(t_data0) * stepT + dT
        s = self._pdfNO_func(t_pdf0) * stepT * self.scale
        tmp_nll += np.sum(s)
        nll = tmp_nll
        return nll
        

    def calc_binnedNLL_IO(self, data, dT) -> float:
        nll = 0
        stepT = 0.01
        # bins with counts
        t_data = np.where(data) 
        n = data[t_data]
        t_data = np.array(t_data) * stepT + self.fitTmin
        t_pdf = t_data + dT 
        s = self._pdfIO_func(t_pdf) * stepT * self.scale
        tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
        tmp_nll = np.sum(tmp_nll) 
        # bins without counts
        t_data0 = np.where(data == 0)
        t_pdf0 = self.fitTmin + np.array(t_data0) * stepT + dT
        s = self._pdfIO_func(t_pdf0) * stepT * self.scale
        tmp_nll += np.sum(s)
        nll = tmp_nll
        return nll

    def calc_Asimov_NLL_NO(self, dT, ty) -> float:
        nll = 0

        stepT = self.Tbinwidth
        binlow, binhig = int(self.fitTmin / stepT) - 1, int(self.fitTmax / stepT) + 1
        for ibin in range(binlow, binhig, 1):
            t_data = stepT * (ibin + 0.5)
            if self.MH == "NO":
                if ty == "MO":
                    n = self._pdfNO_func(t_data) * stepT + self.c14rate * stepT
                elif ty == "Mass":
                    n = self._pdfNO_func0(t_data) * stepT + self.c14rate * stepT
            else:
                if ty == "MO":
                    n = self._pdfIO_func(t_data) * stepT + self.c14rate * stepT
                elif ty == "Mass":
                    n = self._pdfIO_func0(t_data) * stepT + self.c14rate * stepT
            t_pdf = t_data + dT
            s = self._pdfNO_func(t_pdf) * stepT + self.c14rate * stepT
            
            n = n * self.scale
            s = s * self.scale

            if s != 0 and n!=0:
                tmp_nll = s - n + n * np.log(n/s)
                #tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                nll += tmp_nll
            if s != 0 and n == 0:
                tmp_nll = s
                nll += tmp_nll

        return nll



    def calc_Asimov_NLL_NO2D(self, dT, ty) -> float:
        nll = 0
        stepT = self.Tbinwidth
        stepE = self.Ebinwidth
        Tbinlow, Tbinhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        Ebinlow, Ebinhig = int(self.Ethr / stepE), int(self.fitEmax / stepE)
        for it in tqdm(range(Tbinlow, Tbinhig, 1)):
            for iE in range(Ebinlow, Ebinhig, 1):
                t_data = stepT * (it + 0.5)
                E_data = stepE * (iE + 0.5)

                if self.MH == "NO":
                    if ty == "MO":
                        n = self._pdf2DNO_func(t_data, E_data) * stepT * stepE
                    elif ty == "Mass":
                        n = self._pdf2DNO_func0(t_data, E_data) * stepT * stepE
                else:
                    if ty == "MO":
                        n = self._pdf2DIO_func(t_data, E_data) * stepT * stepE
                    elif ty == "Mass":
                        n = self._pdf2DIO_func0(t_data, E_data) * stepT * stepE


                t_pdf = t_data + dT
                s = self._pdf2DNO_func(t_pdf, E_data) * stepT * stepE

                n = n * self.scale
                s = s * self.scale

                if s != 0 and n!=0:
                    #tmp_nll = s - n + n * np.log(n/s)
                    tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                    nll += tmp_nll
                if s != 0 and n == 0:
                    tmp_nll = s
                    nll += tmp_nll

        return -nll
            

    def calc_Asimov_NLL_IO(self, dT, ty) -> float:
        nll = 0

        stepT = self.Tbinwidth
        binlow, binhig = int(self.fitTmin / stepT) - 1, int(self.fitTmax / stepT) + 1
        for ibin in range(binlow, binhig, 1):
            t_data = stepT * (ibin + 0.5)
            if self.MH == "NO":
                if ty == "MO":
                    n = self._pdfNO_func(t_data) * stepT + self.c14rate * stepT
                elif ty == "Mass":
                    n = self._pdfNO_func0(t_data) * stepT + self.c14rate * stepT
            else:
                if ty == "MO":
                    n = self._pdfIO_func(t_data) * stepT + self.c14rate * stepT
                elif ty == "Mass":
                    n = self._pdfIO_func0(t_data) * stepT + self.c14rate * stepT
            t_pdf = t_data + dT
            s = self._pdfIO_func(t_pdf) * stepT + self.c14rate * stepT
            
            n = n * self.scale
            s = s * self.scale

            if s != 0 and n!=0:
                tmp_nll = s - n + n * np.log(n/s)
                #tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                nll += tmp_nll
            if s != 0 and n == 0:
                tmp_nll = s
                nll += tmp_nll

        return nll


    def calc_Asimov_NLL_IO2D(self, dT, ty) -> float:
        nll = 0
        stepT = self.Tbinwidth
        stepE = self.Ebinwidth
        Tbinlow, Tbinhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        Ebinlow, Ebinhig = int(self.Ethr / stepE), int(self.fitEmax / stepE)
        for it in tqdm(range(Tbinlow, Tbinhig, 1)):
            for iE in range(Ebinlow, Ebinhig, 1):
                t_data = stepT * (it + 0.5)
                E_data = stepE * (iE + 0.5)

                if self.MH == "NO":
                    if ty == "MO":
                        n = self._pdf2DNO_func(t_data, E_data) * stepT * stepE
                    elif ty == "Mass":
                        n = self._pdf2DNO_func0(t_data, E_data) * stepT * stepE
                else:
                    if ty == "MO":
                        n = self._pdf2DIO_func(t_data, E_data) * stepT * stepE
                    elif ty == "Mass":
                        n = self._pdf2DIO_func0(t_data, E_data) * stepT * stepE

                t_pdf = t_data + dT
                s = self._pdf2DIO_func(t_pdf, E_data) * stepT * stepE

                n = n * self.scale
                s = s * self.scale

                if s != 0 and n!=0:
                    #tmp_nll = s - n + n * np.log(n/s)
                    tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                    nll += tmp_nll
                if s != 0 and n == 0:
                    tmp_nll = s
                    nll += tmp_nll

        return -nll
            





    def calc_NLL_NO_withbkg(self, data, dT, Tbkg) -> float:
        """
        calculate NLL for given dataset with background and PDF, where PDF is shifted by dT.
        """
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        alldata = []
        for i in data:
            alldata.append(i)
        for i in Tbkg:
            alldata.append(i)

        for i in alldata:
            tmp_nll = np.interp(i+dT, self.pdfNOx, self.pdfNOy) * self.scale + self.c14rate * 1e-3
            if tmp_nll <= 0:
                tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        #bin1 = self.pdfNO.GetXaxis().FindBin(tmin)
        #bin2 = self.pdfNO.GetXaxis().FindBin(tmax)
        #intg = self.pdfNO.Integral(bin1, bin2, "width")

        intg = integrate.quad(self._pdfNO_func, tmin, tmax)[0] * self.scale + (self.fitTmax - self.fitTmin) * self.c14rate * 1e-3
        nll -= intg    
        return -nll

        
    def calc_NLL_IO_withbkg(self, data, dT, level, Tbkg) -> float:
        """
        calculate NLL for given dataset and PDF, where PDF is shifted by dT.
        """
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        alldata = []
        for i in data:
            alldata.append(i)
        for i in Tbkg:
            alldata.append(i)

        for i in alldata:
            tmp_nll = np.interp(i+dT, self.pdfNOx, self.pdfNOy) * self.scale + self.c14rate * 1e-3
            if tmp_nll <= 0:
                tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        #bin1 = self.pdfNO.GetXaxis().FindBin(tmin)
        #bin2 = self.pdfNO.GetXaxis().FindBin(tmax)
        #intg = self.pdfNO.Integral(bin1, bin2, "width")

        intg = integrate.quad(self._pdfNO_func, tmin, tmax)[0] * self.scale + (self.fitTmax - self.fitTmin) * self.c14rate * 1e-3
        nll -= intg    
        return -nll
        




