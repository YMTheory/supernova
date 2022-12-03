from SNNuloader import SNNuloader
from IBD_XS import  IBD_XS
from NuE_XS import NuE_XS
from CEvNS_XS import CEvNS_XS
from toyDetector import  toyDetector
from visible_spectrum import visible_spectrum

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    #loader = SNNuloader(25, "Accretion", "LS220", 10)
    #loader.load()
    det = toyDetector(20000, 0.12, 0.88)
    Na = 6.022e23
    gtoeV = 5.61e32
    mC12 = 12 / Na * gtoeV
    xs = CEvNS_XS(6, 6, mC12)
    dist = 10
    model0 = SNNuloader(25, "Accretion", "LS220", dist)

    chaname = "CEvNS"
    spec = visible_spectrum(chaname, model0, xs, det)

    val = spec.getVisibleEventAtT(0.02, 0.1, 0.105, "NO")
    print("event number per ms : ", val)

