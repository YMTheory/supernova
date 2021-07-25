import numpy as np
from scipy.special import gamma

class EvEneSpec(object):

    m_totE = 5e52   # erg unit
    m_aveE = 12     # MeV unit
    m_dist = 10     # kpc
    gam    = 3 
    m_epsilon = 1

    def __init__(self, totE, aveE, d):
        self.m_totE = totE
        self.m_aveE = aveE
        self.m_dist = d


    def KRJ_dFdE(self, Ev):
        return 3.5e13 / (4*np.pi*self.m_dist**2) * 1 / (self.m_aveE) * Ev**self.gam / gamma(1+self.gam) * ((1+self.gam)/self.m_aveE)**(1+self.gam) * np.exp(-(1+self.gam)*Ev/self.m_aveE)
