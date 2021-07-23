import numpy as np

class detector(object):
   
    mLS = 20   #kton
    fracP = 12e-2
    fracC = 88e-2
    NA = 6.022e23
    molarH = 1   # g/mol
    molarC = 12
    fNumP = mLS*1e9*fracP*NA/molarH
    fNumC = mLS*1e9*fracC*NA/molarC
    fEthr = 0.2   # MeV

    def __init__(self):
        pass

    def getNumberOfProton(self):
        return self.fNumP

    def getNumberOfCarbon(self):
        return self.fNumC

    def getThresholdE(self):
        return self.fEthr

