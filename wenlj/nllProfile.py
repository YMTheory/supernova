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
        return groupNum, deltaT, nllVal, nsig

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()

    if num_datasets > len(lines):
        print(filename)
        print('Error! num_datasets is larger than the number of rows in txt file!')
        return groupNum, deltaT, nllVal, nsig

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
    dist = 3.0
    Ethr = 0.1
    MH = 2
    group = 1
    mNu, minNLL = [], []

    prefix = './dataset/fineSpec/%2.1fkpc/TEvisDATA_'%(dist)
    prefix += 'mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname],MH)

    # Get the number of tests in one group
    resFilename = prefix+'_data%2.1feV_pdf0.0eV_%2.1fkpc_Ethr%2.1fMeV_group%d'%(nuMass1, dist, Ethr, group)
    #resFilename = prefix+'_data%2.1feV_pdf0.0eV_%2.1fkpc_group%d'%(nuMass1, dist, group)
    Tmax = 0.02
    postfix = '_Tmax%3.2f'%Tmax
    resFilename = resFilename + '_1DFitSummary' + postfix + '.txt'

    if os.path.exists(resFilename) == False:
        sys.exit(1)

    file = open(resFilename, 'r')
    lines = file.readlines()
    nStatsPerGroup = len(lines)
    file.close()
    #########################
    #nStatsPerGroup = 100 ## for 2 kpc
    nStatsPerGroup = 500 ## for 2 kpc

    num_datasets = 500
    numass_hypotheses = 19
    xvals = np.linspace(0, 2.0, numass_hypotheses)
    lambdas = np.zeros((num_datasets, numass_hypotheses)) # Empty array which will contain the computed test statistics
    #print(xvals)
    bestfitT = np.zeros((num_datasets, numass_hypotheses))
    nGroups = int(num_datasets/nStatsPerGroup)
    print('nGroups, nStatsPerGroup', nGroups, nStatsPerGroup)

    fig, axs   = plt.subplots(3, 7, figsize=(30, 5))
    #fig2, axs2 = plt.subplots(3, 7, figsize=(30, 5))

    for iPDF in range(numass_hypotheses):
        # special for 2.0 kpc
        #if iPDF==11:
        #    continue
        if iPDF<4:
            continue

        groupNum, deltaT, nllVal, nsig = [], [], [], []
        for iGroup in range(group, nGroups+group):
            nuMass2 = iPDF * 0.1 # eV
            resFilename = prefix+'_data%2.1feV_pdf%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d'%(nuMass1, nuMass2, dist, Ethr, iGroup)
            #resFilename = prefix+'_data%2.1feV_pdf%2.1feV_%2.1fkpc_group%d'%(nuMass1, nuMass2, dist, iGroup)
            resFilename = resFilename + '_1DFitSummary' + postfix + '.txt'

            t1, t2, t3, t4 = getNLL(resFilename, nStatsPerGroup)
            groupNum = groupNum + t1
            deltaT = deltaT + t2
            nllVal = nllVal + t3
            nsig   = nsig   + t4

        for iData in range(num_datasets):
            lambdas[iData, iPDF] = nllVal[iData] * 2 # chi2 ~ 2* (NLL(x) - NLL_loc_min)
            # bestfitT[iData, iPDF] = deltaT[iData]
        
        maxT = np.max(deltaT)
        minT = np.min(deltaT)
        tBins = int( (maxT - minT)/0.0001 )
        print('iPDF, tBins, minT, maxT', iPDF, tBins, minT, maxT)
        if tBins==0:
            continue
        ix = int(iPDF/7)
        iy = iPDF - 7*ix
        axs[ix, iy] = plt.subplot(3, 7, 1+iPDF)
        axs[ix, iy].hist(deltaT, bins=tBins, range=(minT, maxT))
        axs[ix, iy].set_ylabel('at %2.1f eV'%(nuMass2))
        #axs[ix, iy].set_xlabel('$\delta_{T}$ (s)')
        #axs[ix, iy].set_title('%2.1f eV'%(nuMass2))
        # axs[ix, iy] = plt.subplot(3, 7, 1+iPDF)
        # axs[ix, iy].hist(nsig, 20)
        # axs[ix, iy].set_ylabel('at %2.1f eV'%(nuMass2))
    #plt.show()
    picName = './spectra/fineSpec/dT_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_MH%d'%\
              (nuMass1, dist, Ethr, group, MH)
    picName = picName + postfix
    plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')
    #plt.savefig('./spectra/dT_data%2.1feV_%2.1fkpc_group%d.png'%(nuMass1, dist, group), dpi=400, bbox_inches='tight')
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
            elif crossings[iData]>1.5:
                crossings_mask[iData] = False

    plt.plot(xvals,np.ones(len(xvals))*2.706,'--k')
    plt.text(0.5,2.9,'90% confidence threshold (Wilks\')',fontsize=14)
    plt.axis([0.,2.,0.,10.])
    #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
    plt.ylabel('Test statistic')
    plt.xlabel(r'$m_{\nu}$ (eV)')
    picName = './spectra/fineSpec/numass_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_MH%d'%\
              (nuMass1, dist, Ethr, group, MH)
    picName = picName + postfix
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')
    #plt.savefig('./spectra/numass_data%2.1feV_%2.1fkpc_group%d.png'%(nuMass1, dist, group), dpi=400, bbox_inches='tight')

    meanSens = np.mean(crossings[crossings_mask])
    medianSens = np.median(crossings[crossings_mask])
    print(meanSens,'---mean----')
    print(medianSens,'---median----')

    fig3, ax3 = plt.subplots(1,1)
    nCounts, bins, patches = plt.hist(crossings[crossings>0.],bins=np.linspace(0,2,41))
    nCounts = np.array(nCounts)
    ax3.set_ylabel('Counts ({} trials total)'.format(num_datasets))
    ax3.set_xlabel('neutrino mass (eV) at 90% confidence limit')
    plt.plot([meanSens, meanSens], [0., np.max(nCounts)],'--k')
    plt.plot([medianSens, medianSens], [0., np.max(nCounts)],':k')
    #plt.savefig('./spectra/fineSpec/sens_%2.1fkpc_group%d.png'%(dist, group),dpi=400,bbox_inches='tight')
    picName = './spectra/fineSpec/sens_%2.1fkpc_group%d_MH%d'%(dist, group, MH)
    picName = picName + postfix
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')

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
