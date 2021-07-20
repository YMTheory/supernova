import ROOT
import math
import string
import matplotlib.pyplot as plt
import sys, os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

def getValue(strline):
    values = strline.split()
    x1 = int(values[0])
    x2 = float(values[1])
    x3 = float(values[2])
    x4 = float(values[3])
    return (x1, x2, x3, x4)

def func(x, a0, t0, a1):
    return a0*np.power(x-t0, 2) + a1

def getNLL(filename, num_datasets):
    groupNum, deltaT, nllVal, nsig = [], [], [], []
    
    if os.path.exists(filename) == False:
        return False

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()

    if num_datasets > len(lines):
        print('Error! num_datasets is larger than the number of rows in txt file!')
        return False

    for k in range(num_datasets):
        (num, dT, val, nEvt) = getValue( lines[k] )
        groupNum.append(num)
        deltaT.append(dT)
        nllVal.append(val)
        nsig.append(nEvt)
    return groupNum, deltaT, nllVal, nsig
    # deltaT = np.array(deltaT)
    # nllVal = np.array(nllVal)
    # print(deltaT)
    # print(nllVal)
    # locMinNLL = np.min(nllVal)

if __name__ == "__main__":
    
    nuMass1 = 0.0
    #nuMass2 = 0.2
    modelNum = 82503
    chaname = 1
    chaName = ['nup','nue','IBD']
    dist = 10.0
    Ethr = 0.1
    MH = 0
    group = 1
    mNu, minNLL = [], []

    prefix = "./dataset/%2.1fkpc/TEvisDATA_"%(dist)
    prefix += 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)

    num_datasets = 400
    numass_hypotheses = 21
    xvals = np.linspace(0, 2.0, numass_hypotheses)
    lambdas = np.zeros((num_datasets, numass_hypotheses)) # Empty array which will contain the computed test statistics
    print(xvals)

    bestfitT = np.zeros((num_datasets, numass_hypotheses))
    
    #plt.figure()
    fig, axs   = plt.subplots(3, 7, figsize=(30, 5))
    #fig2, axs2 = plt.subplots(3, 7, figsize=(30, 5))

    for iPDF in range(numass_hypotheses):
        nuMass2 = iPDF * 0.1 # eV
        resFilename = prefix+'_data%2.1feV_pdf%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d'%(nuMass1, nuMass2, dist, Ethr, group)
        resFilename = resFilename + '_1DFitSummary.txt'

        groupNum, deltaT, nllVal, nsig = getNLL(resFilename, num_datasets)

        for iData in range(num_datasets):
            lambdas[iData, iPDF] = nllVal[iData]
            # bestfitT[iData, iPDF] = deltaT[iData]
        
        ix = int(iPDF/7)
        iy = iPDF - 7*ix
        axs[ix, iy] = plt.subplot(3, 7, 1+iPDF)
        axs[ix, iy].hist(deltaT, 20)
        axs[ix, iy].set_ylabel('at %2.1f eV'%(nuMass2))
        #axs[ix, iy].set_xlabel('$\delta_{T}$ (s)')
        #axs[ix, iy].set_title('%2.1f eV'%(nuMass2))
        # axs[ix, iy] = plt.subplot(3, 7, 1+iPDF)
        # axs[ix, iy].hist(nsig, 20)
        # axs[ix, iy].set_ylabel('at %2.1f eV'%(nuMass2))
    #plt.show()
    plt.savefig('./spectra/dT_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d.png'%(nuMass1, dist, Ethr, group), dpi=400, bbox_inches='tight')
    # plt.savefig('./dataset/nsig_data%2.1feV_%2.1fkpc_group%d.png'%(nuMass1, dist, group))

    # plot the nll profile
    crossings = np.zeros(num_datasets) - 1.
    crossings_mask = np.ones(num_datasets,dtype=bool)
    fig2, axs2 = plt.subplots(figsize=(6, 4))
    for iData in range(num_datasets):
        locMin = np.min(lambdas[iData])
        nllShifted = lambdas[iData] - locMin
        plt.plot(xvals,nllShifted,'-',label='ToyData{}'.format(iData),linewidth=0.5)

        if iData%20==0:
            print(iData)
        # only fit around the threhsold...
        mask = (nllShifted>=0.)&(nllShifted<7.)
        #print(mask)

        if len(xvals[mask]) > 0:
            p = np.polyfit(xvals[mask],nllShifted[mask],2.)
            crossings[iData] = (-p[1] + np.sqrt( p[1]**2 - 4*(p[0])*(p[2]-2.706) ))/(2*p[0])
            if math.isnan(crossings[iData])==True:
                crossings_mask[iData] = False

    print(crossings)
    plt.plot(xvals,np.ones(len(xvals))*2.706,'--k')
    plt.text(0.5,2.9,'90% confidence threshold (Wilks\')',fontsize=14)
    plt.axis([0.,2.,0.,10.])
    #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
    plt.ylabel('Test statistic')
    plt.xlabel(r'$m_{\nu}$ (eV)')
    plt.savefig('./spectra/numass_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d.png'%(nuMass1, dist, Ethr, group), dpi=400, bbox_inches='tight')

    print(np.mean(crossings[crossings_mask]),'---mean----')
    print(np.median(crossings[crossings_mask]),'---median----')

    fig3, ax3 = plt.subplots(1,1)
    hteststats = plt.hist(crossings[crossings>0.],bins=np.linspace(0,2,40))
    ax3.set_ylabel('Counts ({} trials total)'.format(num_datasets))
    ax3.set_xlabel('neutrino mass (eV) at 90% confidence limit')
    plt.savefig('./spectra/sens_%2.1fkpc.png'%dist,dpi=400,bbox_inches='tight')

    #fig3, ax3 = plt.subplots(1,1)

    # plot the bestfit dT
    # plt.figure()
    # for iPDF in range(numass_hypotheses):
    #     nuMass2 = iPDF * 0.1    # eV
    #     #plotId = '25%d'%iPDF
    #     plt.subplot(2, 5, 1+iPDF)
    #     plt.hist(deltaT)
    #     plt.yscale('Num. of datasets')
    #     plt.xscale('$\delta_{T}$ (s)')
    #     plt.title('%2.1f eV'%(nuMass2))
