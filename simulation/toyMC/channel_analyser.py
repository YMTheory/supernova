import numpy as np
import os, sys
import uproot as up
from collections import Counter
import scipy.integrate as integrate
from scipy.special import gamma
#import ROOT
import time

class channel :
    """
    Define a new class describing a specific detection channel in JUNO.
    """
    
    def __init__(self, name:str, MH:str, model:str, modelNo:int, Ethr:float, dist=10, fitTmin=10, fitTmax=50, fileNo=0, exp="") -> None:
        self.name   = name
        self.MH     = MH
        self.model  = model
        self.modelNo = modelNo
        self.Ethr   = Ethr
        self.dist   = dist
        self.scale  = 10. * 10. / dist / dist

        self.fitTmin = fitTmin
        self.fitTmax = fitTmax

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
        self.pdfNOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_NO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        self.pdfIOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_IO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        #self.pdfNOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs/1D/{model}{modelNo}_PDF_NO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        #self.pdfIOfile = f"/junofs/users/miaoyu/supernova/simulation/C++/PDFs/1D/{model}{modelNo}_PDF_IO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"

        ####### Datasets and PDFs
        self.data_array = None
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

        self.glow       = None
        self.ghig       = None
        self.c14rate    = 0

        self.binwidth   = 0.01

    def setNevtPerFile(self, N) -> None:
        self.NevtPerFile = N

    def getNevtPerFile(self) -> int:
        return self.NevtPerFile
    
    def setDataFilePath(self, datafile) -> None:
        self.datafile = datafile

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

    def _load_pdf(self) -> None:
        """
        Load TH1 PDFs from root files for both NO and IO cases.
        """
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
        
        
    def get_one_event(self, event_id:int):
        """
        get a single event with ID == event_id.
        """
        event_id = event_id - self.startEvt
        return self.data_array[event_id]
    
    def get_one_binned_event(self, event_id:int):
        event_id = event_id - self.startEvt
        return self.binned_data_array[event_id]
    
    
    def _pdfNO_func(self, t):
        return np.interp(t, self.pdfNOx, self.pdfNOy)
    
    def _pdfIO_func(self, t):
        return np.interp(t, self.pdfIOx, self.pdfIOy)

    
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
                tmp_nll = 1e-10
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
                tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        #bin1 = self.pdfIO.GetXaxis().FindBin(tmin)
        #bin2 = self.pdfIO.GetXaxis().FindBin(tmax)
        #intg = self.pdfIO.Integral(bin1, bin2, "width")
        
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


    def calc_Asimov_NLL_NO(self, dT) -> float:
        nll = 0

        stepT = self.binwidth
        binlow, binhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        for ibin in range(binlow, binhig, 1):
            t_data = stepT * ibin
            if self.MH == "NO":
                n = self._pdfNO_func(t_data) * stepT
            else:
                n = self._pdfIO_func(t_data) * stepT

            t_pdf = t_data + dT
            s = self._pdfNO_func(t_pdf) * stepT

            n = n * self.scale
            s = s * self.scale

            if s != 0 and n!=0:
                #tmp_nll = s - n + n * np.log(n/s)
                tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                nll += tmp_nll
            if s != 0 and n == 0:
                tmp_nll = s
                nll += tmp_nll

        return nll

    def calc_Asimov_NLL_IO(self, dT) -> float:
        nll = 0

        stepT = self.binwidth
        binlow, binhig = int(self.fitTmin / stepT), int(self.fitTmax / stepT)
        for ibin in range(binlow, binhig, 1):
            t_data = stepT * ibin
            if self.MH == "NO":
                n = self._pdfNO_func(t_data) * stepT
            else:
                n = self._pdfIO_func(t_data) * stepT

            t_pdf = t_data + dT
            s = self._pdfIO_func(t_pdf) * stepT

            n = n * self.scale
            s = s * self.scale

            if s != 0 and n!=0:
                #tmp_nll = s - n + n * np.log(n/s)
                tmp_nll = s - n * np.log(s) + np.log(gamma(n+1))
                nll += tmp_nll
            if s != 0 and n == 0:
                tmp_nll = s
                nll += tmp_nll

        return nll


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
        




