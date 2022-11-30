import numpy as np
from scipy.special import gamma

class SNNuloader :

    def __init__(self, mass, phase, EoS, d) -> None:
        self.mass   = mass
        self.phase  = phase
        self.EoS    = EoS
        self.d      = d

        self.Lnue           = None
        self.Lnuebar        = None # foe/s -> 10^51 erg/s, P.S. 1 erg = 6.24151e5 MeV
        self.Lnux           = None
        self.Tnue           = None
        self.Tnuebar        = None # s
        self.Tnux           = None
        self.aveEnue        = None # MeV
        self.aveEnux        = None # MeV
        self.aveEnuebar     = None # MeV
        self.aveE2nue       = None # MeV2
        self.aveE2nux       = None # MeV2
        self.aveE2nuebar    = None # MeV2
        self.alphanue       = None
        self.alphanuebar    = None
        self.alphanux       = None

    def load(self) -> None:
        path = f"/junofs/users/miaoyu/supernova/simulation/SN-data-share/timedata/{self.phase}/{self.EoS}/s{self.mass:.1f}co/timedata/"
        filename = path + "neutrino_signal_nu_e"
        print(f"filename for nue : {filename}")
        arr = np.loadtxt(filename)
        self.Tnue = arr[:, 0]
        self.Lnue = arr[:, 1]
        self.aveEnue = arr[:, 2]
        self.aveE2nue = arr[:, 3]
        filename = path + "neutrino_signal_nu_x"
        print(f"filename for nux : {filename}")
        arr = np.loadtxt(filename)
        self.Tnux = arr[:, 0]
        self.Lnux = arr[:, 1]
        self.aveEnux = arr[:, 2]
        self.aveE2nux = arr[:, 3]
        filename = path + "neutrino_signal_nubar_e"
        print(f"filename for nuebar : {filename}")
        arr = np.loadtxt(filename)
        self.Tnuebar = arr[:, 0]
        self.Lnuebar = arr[:, 1]
        self.aveEnuebar = arr[:, 2]
        self.aveE2nuebar = arr[:, 3]

        self.alphanue = 1./(self.aveE2nue/np.power(self.aveEnue, 2)-1) - 1
        self.alphanuebar = 1./(self.aveE2nuebar/np.power(self.aveEnuebar, 2)-1) - 1
        self.alphanux = 1./(self.aveE2nux/np.power(self.aveEnux, 2)-1) - 1

    def reload(self)->None:
        if self.EoS == "Shen":
            eosid = 9
        elif self.EoS == "LS220":
            eosid = 8
        modNo = int(eosid*1e4 + self.mass*100 + 3)
        path = f"/junofs/users/miaoyu/supernova/wenlj/simulation/data/Garching/{modNo}/timedata/"
        filename = path + "neutrino_signal_nu_e"
        print(f"filename for nue : {filename}")
        arr = np.loadtxt(filename)
        self.Tnue = arr[:, 0]
        self.Lnue = arr[:, 1]
        self.aveEnue = arr[:, 2]
        self.aveE2nue = arr[:, 3]
        filename = path + "neutrino_signal_nu_x"
        print(f"filename for nux : {filename}")
        arr = np.loadtxt(filename)
        self.Tnux = arr[:, 0]
        self.Lnux = arr[:, 1]
        self.aveEnux = arr[:, 2]
        self.aveE2nux = arr[:, 3]
        filename = path + "neutrino_signal_nubar_e"
        print(f"filename for nuebar : {filename}")
        arr = np.loadtxt(filename)
        self.Tnuebar = arr[:, 0]
        self.Lnuebar = arr[:, 1]
        self.aveEnuebar = arr[:, 2]
        self.aveE2nuebar = arr[:, 3]

        self.alphanue = 1./(self.aveE2nue/np.power(self.aveEnue, 2)-1) - 1
        self.alphanuebar = 1./(self.aveE2nuebar/np.power(self.aveEnuebar, 2)-1) - 1
        self.alphanux = 1./(self.aveE2nux/np.power(self.aveEnux, 2)-1) - 1

    def getLumAtTime(self, t, flavorID):
        """
        input t: unit sec
        return neutrino luminosity at time t: unit erg
        """
        if flavorID == 0: # nu_e
            if t < self.Tnue[0] or t > self.Tnue[-1]:
                return 0
            else:
                return np.interp(t, self.Tnue, self.Lnue)
        if flavorID == 1: # nu_e
            if t < self.Tnuebar[0] or t > self.Tnuebar[-1]:
                return 0
            else:
                return np.interp(t, self.Tnuebar, self.Lnuebar)
        if flavorID in [2, 3, 4, 5]: # nu_e
            if t < self.Tnux[0] or t > self.Tnux[-1]:
                return 0
            else:
                return np.interp(t, self.Tnux, self.Lnux)

    def getAverageEAtTime(self, t, flavorID):
        """
        input t: unit sec
        return avergaeE of nu at time t: unit MeV
        """
        if flavorID == 0: # nu_e
            if t < self.Tnue[0] or t > self.Tnue[-1]:
                return 0
            else:
                return np.interp(t, self.Tnue, self.aveEnue)
        if flavorID == 1: # nu_e
            if t < self.Tnuebar[0] or t > self.Tnuebar[-1]:
                return 0
            else:
                return np.interp(t, self.Tnuebar, self.aveEnuebar)
        if flavorID in [2, 3, 4, 5]: # nu_e
            if t < self.Tnux[0] or t > self.Tnux[-1]:
                return 0
            else:
                return np.interp(t, self.Tnux, self.aveEnux)

    def getAlphaAtTime(self, t, flavorID):
        """
        input t: unit sec
        return avergaeE of nu at time t: unit MeV
        """
        if flavorID == 0: # nu_e
            if t < self.Tnue[0] or t > self.Tnue[-1]:
                return 0
            else:
                return np.interp(t, self.Tnue, self.alphanue)
        if flavorID == 1: # nu_e
            if t < self.Tnuebar[0] or t > self.Tnuebar[-1]:
                return 0
            else:
                return np.interp(t, self.Tnuebar, self.alphanuebar)
        if flavorID in [2, 3, 4, 5]: # nu_e
            if t < self.Tnux[0] or t > self.Tnux[-1]:
                return 0
            else:
                return np.interp(t, self.Tnux, self.alphanux)

    def getEventAtTime(self, t, E, flavorID):
        """
        input t: unit sec
        input E: unit MeV
        return event number of nu at time t: unit s-1*MeV-1
        """
        lum = self.getLumAtTime(t, flavorID)
        aveE = self.getAverageEAtTime(t, flavorID)
        if aveE == 0:
            return 0
        alpha = self.getAlphaAtTime(t, flavorID)
        indexG = 6.24151e56
        flux = indexG * (lum/aveE) * (np.power(E, alpha) / gamma(1+alpha)) * np.power((alpha+1)/aveE, alpha+1) * np.exp(-(alpha+1)*E/aveE)
        #print(f"flavor {flavorID}, {t}, {E}, {lum}, {aveE}, {alpha}, {flux/indexG}")
        if flavorID == 0:
            if t < self.Tnue[0] or t > self.Tnue[-1]:
                return 0
            else:
                return flux
        if flavorID == 1:
            if t < self.Tnuebar[0] or t > self.Tnuebar[-1]:
                return 0
            else:
                return flux
        if flavorID in [2, 3, 4 ,5]:
            if t < self.Tnux[0] or t > self.Tnux[-1]:
                return 0
            else:
                return flux

    def getEventAtEarthAtTime(self, t, E, mo, flavorID):
        index = 1./(4*np.pi*np.power(3.086e21, 2)) / self.d**2
        fluence = 0
        if mo == "noosc":
            fluence = self.getEventAtTime(t, E, flavorID)
        else:
            p       = 0 
            pbar    = 0
            if mo == "no":
                p       = 0.022
                pbar    = 0.687
            if mo == "io":
                p       = 0.291
                pbar    = 0.022
            
            if flavorID == 0:
                fluence = p * self.getEventAtTime(t, E, 0) + (1 - p) * self.getEventAtTime(t, E, 2)
            if flavorID == 1:
                fluence = pbar * self.getEventAtTime(t, E, 1) + (1 - pbar) * self.getEventAtTime(t, E, 3)
                #print(pbar, self.getEventAtTime(t, E, 1), self.getEventAtTime(t, E, 3), fluence, index, fluence*index )
            if flavorID in [2, 4]:
                fluence = 0.5 * (1 - p) * self.getEventAtTime(t, E, 0) + 0.5 * p * self.getEventAtTime(t, E, 2)
            if flavorID in [3, 5]:
                fluence = 0.5 * (1 - pbar) * self.getEventAtTime(t, E, 1) + 0.5 * pbar * self.getEventAtTime(t, E, 3)


        return index * fluence

