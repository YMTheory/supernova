import numpy as np
from astropy import units as u

class toyDetector :

    def __init__(self, mass, fracH, fracC)->None:
        self.mass = mass * 1e6
        self.fracH = fracH
        self.fracC = fracC

        Na = 6.022e23
        self.nH = self.mass * fracH * Na / 1
        self.nC = self.mass * fracC * Na / 12.
    
