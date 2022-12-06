import numpy as np
import uproot as up

class toyDetector :

    def __init__(self, mass, fracH, fracC)->None:
        self.mass = mass * 1e6
        self.fracH = fracH
        self.fracC = fracC

        Na = 6.022e23
        self.nH = self.mass * fracH * Na / 1
        self.nC = self.mass * fracC * Na / 12.
    
        f = up.open("/junofs/users/miaoyu/supernova/wenlj/simulation/data/quenching/quenchCurve.root")
        g = f["quenchCurve"]
        self.proton_quench_x = g.values()[0]
        self.proton_quench_y = g.values()[1]
        g = f["quenchCurveInv"]
        self.proton_unquench_x = g.values()[0]
        self.proton_unquench_y = g.values()[1]


    #@profile
    def proton_quenchedE(self, T):
        return np.interp(T, self.proton_quench_x, self.proton_quench_y)


    #@profile
    def proton_unquenchedE(self, T):
        return np.interp(T, self.proton_unquench_x, self.proton_unquench_y)

