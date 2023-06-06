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
import logging

class channel :
    """
    Define a new class describing a specific detection channel in JUNO.
    """
    
    def __init__(self, name:str, MH:str, model:str, modelNo:int, Ethr:float, dist=10, fitTmin=-0.02, fitTmax=0.02, fitEmax=80, fileNo=0, exp="", nuMass=0.0) -> None:
        self.name   = name
        self.MH     = MH
        self.model  = model
        self.modelNo = modelNo
        self.Ethr   = Ethr
        self.dist   = dist
        self.scale  = 10. * 10. / dist / dist
        self.bkgscale = 1.
        self.nuMass = nuMass

        self.fitTmin = fitTmin
        self.fitTmax = fitTmax
        self.fitEmax = fitEmax

        self.fileNo  = fileNo
        self.NevtPerFile = 1000

        ## Initialize PDFs:
        ### bkg
        self.c14file = None
        self.c14pdf2Dx = None
        self.c14pdf2Dy = None
        self.c14pdf2D = None
        self.f2dC14 = None

        ### 1D:
        self.datapdffile1D = None
        self.datapdf1Dx = None
        self.datapdf1Dy = None

        self.pdfNOfile1D = None
        self.pdfNOhist1D = None
        self.pdfNO1Dx     = None
        self.pdfNO1Dbins  = None
        self.pdfNO1Dy     = None
        self.pdfIOfile1D = None
        self.pdfIOhist1D = None
        self.pdfIO1Dx     = None
        self.pdfIO1Dbins  = None
        self.pdfIO1Dy     = None

        ### 2D:
        self.datapdffile2D = None
        self.datapdf2DT    = None
        self.datapdf2DE    = None
        self.datapdf2D    = None
        self.f2ddatapdf   = None

        self.pdfNOfile2D = None
        self.pdf2DNOT     = None
        self.pdf2DNOE     = None
        self.pdf2DNO     = None
        self.f2dNO      = None
        self.pdfIOfile2D = None
        self.pdf2DIOT     = None
        self.pdf2DIOE     = None
        self.pdf2DIO     = None
        self.f2dIO      = None
        
        ## Initialize data
        self.datafile1D = None
        self.datafile2D = None

        self.data_array = None
        self.dataT_array = None
        self.dataE_array = None
        self.nentries   = 0
        self.startEvt   = 0
        self.endEvt     = self.NevtPerFile

        self.load1DData_flag = False
        self.load2DData_flag = False
        self.load1DPDF_flag = False
        self.load2DPDF_flag = False
        self.load1DdataPDF_flag = False
        self.load2DdatdataaPDF_flag = False

        self.c14rate    = 0
        self.c14concentration = 0.0

        self.Tbinwidth   = 0.00001
        self.Ebinwidth  = 0.2
        if self.name == "pES":
            self.Ebinwidth = 0.01


    def setBkgScale(self, s):
        self.bkgscale = s

    def setScale(self, s):
        self.scale = s

    def setMO(self, MO):
        self.MH = MO

    def setC14rate(self, rate):
        self.c14rate = rate

    def setC14concentration(self, c):
        self.c14concentration = c

    def setC14PdfFile2D(self, pdffile) -> None:
        self.c14file = pdffile

    def setDistance(self, dist):
        self.dist = dist

    def setNevtPerFile(self, N) -> None:
        self.NevtPerFile = N

    def getNevtPerFile(self) -> int:
        return self.NevtPerFile
    
    def setDataFile1D(self, datafile) -> None:
        self.datafile1D = datafile

    def setDataFile2D(self, datafile) -> None:
        self.datafile2D = datafile

    def setDataPdfFile1D(self, pdffile) -> None:
        self.datapdffile1D = pdffile

    def setDataPdfFile2D(self, pdffile) -> None:
        self.datapdffile2D = pdffile

    def setNOPdfFile1D(self, pdffile) -> None:
        self.pdfNOfile1D = pdffile

    def setIOPdfFile1D(self, pdffile) -> None:
        self.pdfIOfile1D = pdffile

    def setNOPdfFile2D(self, pdffile) -> None:
        self.pdfNOfile2D = pdffile

    def setIOPdfFile2D(self, pdffile) -> None:
        self.pdfIOfile2D = pdffile

    def setStartEvtId(self, idd):
        self.startEvt = idd
        
    def setEndEvtId(self, idd):
        self.endEvt = idd


    def _load_data1D(self) -> None:
        """
        Load ttree ak array data from root files.
        """
        try:
            with up.open(self.datafile1D) as f:
                print("--------------------------------------------------------------")
                print(f"Load datafile of channel {self.name} from ->")
                print(self.datafile1D)
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
                self.load1DData_flag = True
        except FileExistsError:
            print(f"The data file {self.datafile1D} dose not exist! :(")
            sys.exit(-1)


    def _load_data2D(self) -> None:
        """
        Load ttree ak array data from root files.
        """
        try:
            with up.open(self.datafile2D) as f:
                logging.debug("--------------------------------------------------------------")
                logging.debug(f"Load 2D datafile2D of channel {self.name} from ->")
                logging.debug(self.datafile2D)
                logging.debug("--------------------------------------------------------------")
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
                logging.debug(f"\n*****Total loaded event number = {evtNum} from Event {self.startEvt} to {self.endEvt} in channel {self.name}, where the average signal number in one event is {Nave} and sigma is {Nsigma}.\n")
                self.dataT_array = nuTime
                self.dataE_array = nuEnergy
                self.nentries = evtNum

                self.load2DPDF_flag = True
        except FileExistsError:
            lbEogging.critical(f"The data file {self.datafile2D} dose not exist! :(")
            sys.exit(-1)


    def _load_pdf1D(self) -> None:
        """
        Load TH1 PDFs from root files for both NO and IO cases.
        """
        logging.debug(f"\n ========= Load {self.name} 1D PDF ========= \n")
        try:
            logging.debug(self.pdfNOfile1D)
            f = up.open(self.pdfNOfile1D)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            print(f"The pdf file {self.pdfNOfile1D} dose not exist! :(")
            sys.exit(-1)
            
        self.pdfNO1Dbins = tmp_h1.to_numpy()[1]
        #self.pdfNOhist1D = tmp_h1.to_hist()
        axis = tmp_h1.axis()    
        self.pdfNO1Dx = axis.centers()
        self.pdfNO1Dy = tmp_h1.values()
        
        try:
            logging.debug(self.pdfIOfile1D)
            f = up.open(self.pdfIOfile1D)
            tmp_h1 = f["h1"] 
        except FileNotFoundError:
            print(f"The pdf file {self.pdfIOfile1D} dose not exist! :(")
            sys.exit(-1)
            
        self.pdfIO1Dbins = tmp_h1.to_numpy()[1]
        #self.pdfIOhist1D = tmp_h1.to_hist()
        axis = tmp_h1.axis()    
        self.pdfIO1Dx = axis.centers()
        self.pdfIO1Dy = tmp_h1.values()
        self.load1DPDF_flag = True

    def _load_datapdf1D(self) -> None:
        """
        Load TH1 PDFs from root files for both NO and IO cases.
        """
        logging.debug(f"\n ========= Load {self.name} 1D data PDF ========= \n")
        try:
            logging.debug(self.datapdffile1D)
            f = up.open(self.datapdffile1D)
            tmp_h1 = f["h1"]
            axis = tmp_h1.axis()
            self.datapdf1Dx = axis.centers()
            self.datapdf1Dy = tmp_h1.values()

            self.load1DdataPDF_flag = True
        except FileNotFoundError:
            print(f"The pdf file {self.datapdffile1D} dose not exist! :(")
            sys.exit(-1)

    def _load_C14pdf2D(self) -> None:
        """
        Load C-14 2D pdf from root files.
        """
        try: 
            logging.debug(f"\n ======= Load {self.c14file} 2D c14 PDF ======== \n")
            f = up.open(self.c14file)
            tmp_h1 = f["c14"]
        except:
            logging.debug(f"The c14 2D pdf file {self.c14file} does not exist! :(")
            sys.exit(-1)

        
        c14concentration0 = 1e-17
        ratio = self.c14concentration / c14concentration0
        xaxis = tmp_h1.axis("x")
        yaxis = tmp_h1.axis("y")
        self.c14pdf2Dx = xaxis.centers()
        self.c14pdf2Dy = yaxis.centers()
        self.c14pdf2D  = tmp_h1.values() * ratio * self.bkgscale
        self.f2dC14 = interpolate.interp2d(self.c14pdf2Dx, self.c14pdf2Dy, self.c14pdf2D.T, kind="linear")

    

    def _load_pdf2D(self) -> None:
        logging.debug(f"\n ========= Load {self.name} 2D PDF ========= \n")
        try:
            logging.debug(self.pdfNOfile2D)
            f = up.open(self.pdfNOfile2D)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            logging.debug(f"The pdf file {self.pdfNOfile2D} does not exist! :(")
            sys.exit(-1)
        xaxis = tmp_h1.axis("x")    
        yaxis = tmp_h1.axis("y")
        self.pdf2DNOT = xaxis.centers()
        self.pdf2DNOE = yaxis.centers()
        self.pdf2DNO  = tmp_h1.values()
        self.f2dNO = interpolate.interp2d(self.pdf2DNOT, self.pdf2DNOE, self.pdf2DNO.T, kind="linear")
             
        try:
            logging.debug(self.pdfIOfile2D)
            f = up.open(self.pdfIOfile2D)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            logging.debug(f"The pdf file {self.pdf2DIOfile} does not exist! :(")
            sys.exit(-1)
        xaxis = tmp_h1.axis("x")    
        yaxis = tmp_h1.axis("y")
        self.pdf2DIOT = xaxis.centers()
        self.pdf2DIOE = yaxis.centers()
        self.pdf2DIO  = tmp_h1.values()
        self.f2dIO = interpolate.interp2d(self.pdf2DIOT, self.pdf2DIOE, self.pdf2DIO.T, kind="linear")

        self.load2DPDF_flag = True
    

    def _load_datapdf2D(self) -> None:
        logging.debug(f"\n ========= Load {self.name} 2D data PDF ========= \n")
        try:
            logging.debug(self.datapdffile2D)
            f = up.open(self.datapdffile2D)
            tmp_h1 = f["h1"]
        except FileNotFoundError:
            logging.debug(f"The pdf file {self.datapdffile2D} does not exist! :(")
            sys.exit(-1)
        xaxis = tmp_h1.axis("x")    
        yaxis = tmp_h1.axis("y")
        self.datapdf2DT = xaxis.centers()
        self.datapdf2DE = yaxis.centers()
        self.datapdf2D  = tmp_h1.values()
        self.f2ddatapdf = interpolate.interp2d(self.datapdf2DT, self.datapdf2DE, self.datapdf2D.T, kind="linear")
    
        self.load2DdatdataaPDF_flag = True

        
    def get_one_event(self, event_id:int):
        """
        get a single event with ID == event_id.
        """
        if not self.load1DData_flag:
            self._load_data_ak()
        event_id = event_id - self.startEvt
        return self.data_array[event_id]
    
    def get_one_event1D(self, event_id:int):
        event_id = event_id - self.startEvt
        return self.dataT_array[event_id]
    
    def get_one_event2D(self, event_id:int):
        event_id = event_id - self.startEvt
        return self.dataT_array[event_id], self.dataE_array[event_id]
    
    def _datapdf_func(self, t):
        return np.interp(t, self.datapdf1Dx, self.datapdf1Dy) 
    
    def _pdfNO_func(self, t):
        return np.interp(t, self.pdfNO1Dx, self.pdfNO1Dy) 
    
    def _pdfIO_func(self, t):
        return np.interp(t, self.pdfIO1Dx, self.pdfIO1Dy)

    def _pdf2DNO_func(self, t, E):
        return self.f2dNO(t, E)[0]

    def _pdf2DIO_func(self, t, E):
        return self.f2dIO(t, E)[0] 

    def _datapdf2D_func(self, t, E):
        return self.f2ddatapdf(t, E)[0]

    def _c14pdf2D_func(self, t, E):
        return self.f2dC14(t, E)[0]

    def getNsigCurrentEvent(self, evtid:int) -> int:
        evtid = evtid - self.startEvt
        return len(self.data_array[evtid])

    def calc_NLL_NO(self, data, dT) -> float:
        """
        calculate NLL for given dataset and PDF, where PDF is shifted by dT.
        """
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for i in data:
            #tmp_nll = np.interp(i+dT, self.pdfNO1Dx, self.pdfNO1Dy) * self.scale
            tmp_nll = _pdfNO_func(i+dT) * self.scale 
            if tmp_nll <= 0:
                continue
            nll += np.log(tmp_nll)
        
        intg = integrate.quad(self._pdfNO_func, tmin, tmax)[0] * self.scale
        nll -= intg    
        return -nll


    def calc_NLL_IO(self, data, dT) -> float:
        """
        calculate NLL for given dataset and PDF, where PDF is shifted by dT.
        """
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for i in data:
            tmp_nll = _pdfIO_func(i+dT) * self.scale
            if tmp_nll <= 0:
                continue
            nll += np.log(tmp_nll)
        
        intg = integrate.quad(self._pdfIO_func, tmin, tmax)[0] * self.scale
        
        nll -= intg 

        return -nll


    def calc_NLL_NO2D(self, dataT, dataE, dT) -> float: ## input data time unit: ms
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        flag = False
        for t, E in zip(dataT, dataE):
            tmp_prob = self._pdf2DNO_func((t+dT), E) * self.scale + self._c14pdf2D_func((t+dT), E)
            #if tmp_prob <= 1e-20:
            #    tmp_prob = 1e-20
            if tmp_prob <= 0:
                flag = True
                tmp_prob = 1e-10
            nll += np.log(tmp_prob)
            #print(dT, t, E, tmp_prob)
        
        intg = integrate.quad(self._pdfNO_func, tmin, tmax)[0] * self.scale
        nll -= intg
        #print(f"Result: {dT}, {intg}, {nll}.")
        return -nll, flag


    def calc_NLL_IO2D(self, dataT, dataE, dT) -> float:
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        flag = False
        for t, E in zip(dataT, dataE):
            #tmp_prob = self._pdf2DIO_func((t+dT), E)  * self.scale
            tmp_prob = self._pdf2DIO_func((t+dT), E) * self.scale + self._c14pdf2D_func((t+dT), E)
            #if tmp_prob <= 1e-20:
            #    tmp_prob = 1e-20
            if tmp_prob <= 0:
                flag = True
                tmp_prob = 1e-10
            nll += np.log(tmp_prob)
        
        intg = integrate.quad(self._pdfIO_func, tmin, tmax)[0] * self.scale
        nll -= intg
        return -nll, flag


    def calc_Asimov_NLL_NO(self, dT) -> float:
        nll = 0

        stepT = self.Tbinwidth
        binlow, binhig = int(self.fitTmin / stepT) - 1, int(self.fitTmax / stepT) + 1
        for ibin in range(binlow, binhig, 1):
            t_data = stepT * (ibin + 0.5)
            n = self._datapdf_func(t_data) * stepT * self.scale #+ self.c14rate * stepT * self.bkgscale
            t_pdf = t_data + dT
            s = self._pdfNO_func(t_pdf) * stepT * self.scale #+ self.c14rate * stepT * self.bkgscale
            
            if s != 0 and n!=0:
                tmp_nll = s - n + n * np.log(n/s)
                #tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                nll += tmp_nll
            if s != 0 and n == 0:
                tmp_nll = s
                nll += tmp_nll

        return nll


    def calc_Asimov_NLL_NO2D(self, dT) -> float:
        nll = 0
        test_n = 0
        stepT = self.Tbinwidth
        stepE = self.Ebinwidth
        Tbinlow, Tbinhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        Ebinlow, Ebinhig = int(self.Ethr / stepE), int(self.fitEmax / stepE)
        for it in tqdm(range(Tbinlow, Tbinhig, 1)):
            for iE in range(Ebinlow, Ebinhig, 1):
                t_data = stepT * (it + 0.5)
                E_data = stepE * (iE + 0.5)
                n = self._datapdf2D_func(t_data, E_data) * stepT * stepE
                test_n += n

                t_pdf = t_data + dT
                s = self._pdf2DNO_func(t_pdf, E_data) * stepT * stepE

                n = n * self.scale
                s = s * self.scale

                if s != 0 and n!=0:
                    tmp_nll = s - n + n * np.log(n/s)
                    #tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                    nll += tmp_nll
                if s != 0 and n == 0:
                    tmp_nll = s
                    nll += tmp_nll
        
        print(f"Total measured event number = {test_n}.")
        return nll


    def calc_Asimov_NLL_NO2D_withC14(self, dT) -> float:
        nll = 0
        stepT, stepE = self.Tbinwidth, self.Ebinwidth
        Tbinlow, Tbinhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        Ebinlow, Ebinhig = int(self.Ethr / stepE), int(self.fitEmax / stepE)
        for it in tqdm(range(Tbinlow, Tbinhig, 1)):
            for iE in range(Ebinlow, Ebinhig, 1):
                t_data = stepT * (it + 0.5)
                E_data = stepE * (iE + 0.5)
                n_sig = self._datapdf2D_func(t_data, E_data) * stepT * stepE

                t_pdf = t_data + dT
                s_sig = self._pdf2DNO_func(t_pdf, E_data) * stepT * stepE

                n_sig = n_sig * self.scale
                s_sig = s_sig * self.scale

                ## Add C14 into fitter:
                n_bkg = self._c14pdf2D_func(t_data, E_data) * stepT * stepE
                s_bkg = self._c14pdf2D_func(t_pdf, E_data) * stepT * stepE
            
                s = s_sig + s_bkg
                n = n_sig + n_bkg

                if s != 0 and n!=0:
                    tmp_nll = s - n + n * np.log(n/s)
                    #tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                    nll += tmp_nll
                if s != 0 and n == 0:
                    tmp_nll = s
                    nll += tmp_nll
        
        return nll

            

    def calc_Asimov_NLL_IO(self, dT) -> float:
        nll = 0
        stepT = self.Tbinwidth
        binlow, binhig = int(self.fitTmin / stepT) - 1, int(self.fitTmax / stepT) + 1
        for ibin in range(binlow, binhig, 1):
            t_data = stepT * (ibin + 0.5)
            n = self._datapdf_func(t_data) * stepT * self.scale #+ self.c14rate * stepT * self.bkgscale
            t_pdf = t_data + dT
            s = self._pdfIO_func(t_pdf) * stepT * self.scale #+ self.c14rate * stepT * self.bkgscale
            
            if s != 0 and n!=0:
                tmp_nll = s - n + n * np.log(n/s)
                #tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                nll += tmp_nll
            if s != 0 and n == 0:
                tmp_nll = s
                nll += tmp_nll

        return nll



    def calc_Asimov_NLL_IO2D(self, dT) -> float:
        nll = 0
        stepT = self.Tbinwidth
        stepE = self.Ebinwidth
        Tbinlow, Tbinhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        Ebinlow, Ebinhig = int(self.Ethr / stepE), int(self.fitEmax / stepE)
        for it in tqdm(range(Tbinlow, Tbinhig, 1)):
            for iE in range(Ebinlow, Ebinhig, 1):
                t_data = stepT * (it + 0.5)
                E_data = stepE * (iE + 0.5)
                n = self._datapdf2D_func(t_data, E_data) * stepT * stepE

                t_pdf = t_data + dT
                s = self._pdf2DIO_func(t_pdf, E_data) * stepT * stepE

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
            

    def calc_Asimov_NLL_IO2D_withC14(self, dT) -> float:
        nll = 0
        stepT, stepE = self.Tbinwidth, self.Ebinwidth
        Tbinlow, Tbinhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        Ebinlow, Ebinhig = int(self.Ethr / stepE), int(self.fitEmax / stepE)
        for it in tqdm(range(Tbinlow, Tbinhig, 1)):
            for iE in range(Ebinlow, Ebinhig, 1):
                t_data = stepT * (it + 0.5)
                E_data = stepE * (iE + 0.5)
                n_sig = self._datapdf2D_func(t_data, E_data) * stepT * stepE

                t_pdf = t_data + dT
                s_sig = self._pdf2DIO_func(t_pdf, E_data) * stepT * stepE

                n_sig = n_sig * self.scale
                s_sig = s_sig * self.scale

                ## Add C14 into fitter:
                n_bkg = self._c14pdf2D_func(t_data, E_data) * stepT * stepE
                s_bkg = self._c14pdf2D_func(t_pdf, E_data) * stepT * stepE
            
                s = s_sig + s_bkg
                n = n_sig + n_bkg

                if s != 0 and n!=0:
                    tmp_nll = s - n + n * np.log(n/s)
                    #tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                    nll += tmp_nll
                if s != 0 and n == 0:
                    tmp_nll = s
                    nll += tmp_nll
        
        return nll

