import ROOT
import math
import string
import matplotlib.pyplot as plt
import sys, os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)
    nuMass = 0.

    fig = plt.subplots(1, 2, figsize=(5.4, 3.6))
    plt.subplots_adjust(left=0.12, top= 0.97, right = 0.97, bottom = 0.14, wspace = 0.01, hspace = 0.01)

    axs = plt.subplot(1, 2, 1)
    #axs = set_xmargin(0.)
    dist = [2., 3., 5., 10.]
    sens1 = [0.55, 0.56, 0.60, 0.66]
    plt.plot( dist, sens1, 'ro-', label='Garching 82503')
    plt.axis([1., 11., 0., 1.])
    plt.xlabel('SN distance (kpc)')
    #axs.set_xscale('log')
    plt.ylabel(r'$m_{\nu}$ median sensitivity (eV)')
    plt.legend(loc='upper left',fontsize=10,ncol=2,labelspacing=0.1)

    axs = plt.subplot(1, 2, 2)
    stat = [1., 4., 25., 100.]
    sens2 = [0.66, 0.44, 0.26, 0.20]
    nominalDist = 10.0
    plt.plot( stat, sens2, 'bs-', label='Garching 82503, 10 kpc')
    plt.axis([0.5, 200., 0., 1.])
    plt.xlabel('N times of nominal statistics')
    #plt.yticks([])
    #axs.yaxis.set_ticks_position('right')
    #axs.yaxis.tick_right()
    for label in axs.get_yticklabels():
        label.set_visible(False)
    axs.set_xscale('log')
    #plt.ylabel(r'$m_{\nu}$ Median sensitivity (eV)')
    plt.legend(loc='upper left',fontsize=10,ncol=2,labelspacing=0.1)
    #plt.show()
    plt.savefig('./spectra/fineSpec/sens_data%2.1feV.png'%(nuMass), dpi=400, bbox_inches='tight')
    plt.savefig('./spectra/fineSpec/sens_data%2.1feV.pdf'%(nuMass), dpi=400, bbox_inches='tight')
