import numpy as np
import os, sys
import uproot as up
from collections import Counter
import scipy.integrate as integrate

class channel :
    """
    Define a new class describing a specific detection channel in JUNO.
    """
    
    def __init__(self, name:str, MH:str, model:str, modelNo:int, Ethr:float, dist=10, fitTmin=10, fitTmax=50, fileNo=0) -> None:
        self.name   = name
        self.MH     = MH
        self.model  = model
        self.Ethr   = Ethr
        self.dist   = dist

        self.fitTmin = fitTmin
        self.fitTmax = fitTmax

        self.fileNo  = fileNo
        
        production_path = "/junofs/users/miaoyu/supernova/production/"
        self.datafile   = f"{production_path}/Data/{dist}kpc/{model}/{model}{modelNo}_{name}_data_{MH}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_new_{fileNo}.root"
        self.pdfNOfile  = f"{production_path}/PDFs/{dist}kpc/{model}/{model}{modelNo}_PDF_NO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"
        self.pdfIOfile  = f"{production_path}/PDFs/{dist}kpc/{model}/{model}{modelNo}_PDF_IO_{dist}kpc_{name}_{Ethr:.2f}MeV_newshortPDF.root"


        ####### Datasets and PDFs
        self.data_array = None
        self.nentries   = 0
        #self.pdfNO      = None
        #self.pdfIO      = None
        self.pdfNOx     = None
        self.pdfNOy     = None
        self.pdfIOx     = None
        self.pdfIOy     = None


    def _load_data(self) -> None:
        """
        Load ttree data from root files.
        """
        try:
            with up.open(self.datafile) as f:
                print(self.datafile)
                evtId  = f["tFixedStat"]["evtID"].array()
                nuTime = f["tFixedStat"]["nuTime1D"].array()

                counter = Counter(evtId)
                nuNum = counter[0]
                evtNum = int(len(nuTime) / nuNum)
                nuTime = np.reshape(nuTime, (evtNum, nuNum))
                nuTime = np.array(nuTime)
        except FileNotFoundError:
            print(f"The data file {self.datafile} dose not exist! :(")
            sys.exit(-1)

        self.data_array = nuTime
        self.nentries   = len(nuTime)


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
        
        
    def get_one_event(self, event_id:int) -> None:
        """
        get a single event with ID == event_id.
        """
        return self.data_array[event_id]
    
    def _pdfNO_func(self, t):
        return np.interp(t, self.pdfNOx, self.pdfNOy)
    
    def _pdfIO_func(self, t):
        return np.interp(t, self.pdfIOx, self.pdfIOy)

    
    def calc_NLL_NO(self, data, dT) -> float:
        """
        calculate NLL for given dataset and PDF, where PDF is shifted by dT.
        """
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for i in data:
            #tmp_nll = self.pdfNO.Interpolate(i + dT)
            tmp_nll = np.interp(i+dT, self.pdfNOx, self.pdfNOy)
            if tmp_nll <= 0:
                tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        #bin1 = self.pdfNO.GetXaxis().FindBin(tmin)
        #bin2 = self.pdfNO.GetXaxis().FindBin(tmax)
        #intg = self.pdfNO.Integral(bin1, bin2, "width")

        intg = integrate.quad(self._pdfNO_func, tmin, tmax)[0]
        nll -= intg    
        return -nll

        
    def calc_NLL_IO(self, data, dT) -> float:
        """
        calculate NLL for given dataset and PDF, where PDF is shifted by dT.
        """
        nll = 0
        tmin, tmax = self.fitTmin + dT, self.fitTmax + dT
        for i in data:
            #tmp_nll = self.pdfIO.Interpolate(i + dT)
            tmp_nll = np.interp(i+dT, self.pdfIOx, self.pdfIOy)
            if tmp_nll <= 0:
                tmp_nll = 1e-10
            nll += np.log(tmp_nll)
        
        #bin1 = self.pdfIO.GetXaxis().FindBin(tmin)
        #bin2 = self.pdfIO.GetXaxis().FindBin(tmax)
        #intg = self.pdfIO.Integral(bin1, bin2, "width")
        
        intg = integrate.quad(self._pdfIO_func, tmin, tmax)[0]
        
        nll -= intg 

        return -nll
        
        
