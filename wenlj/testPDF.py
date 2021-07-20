#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: Liangjian Wen ---> wenlj@ihep.ac.cn
# Created Time : Sat Oct 31
# Generate different (time, Ev) PDFs
"""

import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import TH1D, TH2D, TCanvas
import scipy.stats as st
from scipy.stats import rv_continuous
from array import array
import sys
import os
libSNsimDir = os.getenv('SNCODEDIR') + '/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)
print('Load', libSNsimDir)

if __name__ == "__main__":

    tmin, tmax = -0.1, 0.1
    t1_thr, t2_thr = -0.025, 0.05
    step_t0, step_t1, step_t2 = 0.00005, 0.00005, 0.0005
    npt0 = round( (t1_thr-tmin)/step_t0 )
    npt1 = round( (t2_thr-t1_thr)/step_t1 )
    npt2 = round( (tmax-t2_thr)/step_t2 )
    nbin_time = npt0 + npt1 + npt2
    
    # time binning
    binning_time = array('d', [])
    for ip in range(npt0):
        binning_time.append( tmin+step_t0*ip )
    
    for ip in range(npt0, npt0 + npt1):
        binning_time.append( t1_thr+step_t1*(ip-npt0) )

    for ip in range(npt0 + npt1, npt0 + npt1 + npt2):
        binning_time.append( t2_thr+step_t2*(ip-npt0-npt1) )
    
    binning_time.append(tmax)
    # for debug
    print('npt0, start, end: ', npt0, binning_time[0], binning_time[npt0])
    print('npt1, start, end: ', npt1, binning_time[npt0], binning_time[npt0+npt1])
    print('npt2, start, end: ', npt2, binning_time[npt0+npt1], binning_time[npt0+npt1+npt2])
    print('len(binning_time), nbin_time, last bin:', len(binning_time), nbin_time, binning_time[npt0+npt1+npt2])

