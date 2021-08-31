import os, sys
import numpy as np
import matplotlib.pyplot as plt
import math
import namespace as ns
import numpy.ma as ma

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
        print("%s does not exist !" %filename)
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




if __name__ == "__main__":

    modelNum = int( sys.argv[1] )

    chaName = ['nup','nue','IBD']
    #chaname = snDet.IBD
    chaname = int( sys.argv[2] )
    print("channelName: ", chaname, chaName[chaname])

    nuType = int( sys.argv[3] )
    # -1 for all types; 0 nu_e; 1 anti_nue; 2 nu_x
    nuTypeName = ['nu_e', 'anti_nue', 'nu_x']
    print("neutrino type: ", nuType, nuTypeName[nuType])

    nuMass1 = float( sys.argv[4] )# unit: eV
    print('nuMass for producing dataset : ', nuMass1, ' eV')
    nuMassName1 = int(nuMass1 * 10)

    MHName = ['NoOsc', 'NH', 'IH']
    dataMH = int( sys.argv[5] )
    print('NMO for dataset : ', dataMH, MHName[dataMH])

    fitMH = int( sys.argv[6] )
    print('NMO for fitting : ', fitMH, MHName[fitMH])

    dist = float( sys.argv[7] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[8] )# unit: kpc
    print('fit Ethr: ', Ethr, ' MeV')

    group = int( sys.argv[9] )# group number
    print('group number: ', group)




    ## -----> Load Fit Summary Results <----- ##
    #prefix  = "./dataset/%2.1fkpc/TEvisDATA_" %(dist)
    #prefix += "mod%d_cha%d%s_mh%d" %(modelNum, chaname, chaName[chaname], dataMH) 
    #resFilename = prefix + "_data%2.1feV_mh%d_pdf" %(nuMass1, fitMH)
    #suffix = "eV_%2.1fkpc_Ethr%2.1fMeV_group%d_1DFitSummary_Tmax0.02.txt" %(dist, Ethr, group)
    #print(resFilename+suffix)

    num_datasets = 100 #500
    numass_hypotheses = 30
    nStatsPerGroup = 100 #500
    xvals = np.linspace(0, 3.0, numass_hypotheses)
    lambdas = np.zeros((num_datasets, numass_hypotheses))
    bestfitT = np.zeros((num_datasets, numass_hypotheses))

    
    fig, axs   = plt.subplots(5, 7, figsize=(30, 5))

    for iPDF in range(numass_hypotheses):
        
        groupNum, deltaT, nllVal, nsig = [], [], [], []
        nuMass2 = iPDF * 0.1 #eV
        filename = ns.fitResFileName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)[1]
        #filename = ns.fit2DResFileName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)[1]
        #filename = "./dataset/10.0kpc/TEvisDATA_mod82503_cha1nue_dataMH1_datamNu0.0eV_fitMH1_fitmNu0.0eV_10.0kpc_Ethr0.1MeV_group880_1DFitSummary.txt"
        print(filename)
        if not os.path.exists(filename) :
            print("%s does not exist!")
            continue
        t1, t2, t3, t4 = getNLL(filename, nStatsPerGroup)
        groupNum = groupNum + t1
        deltaT = deltaT + t2
        nllVal = nllVal + t3
        nsig   = nsig   + t4

        for iData in range(num_datasets):
            lambdas[iData, iPDF] = nllVal[iData] * 2

        maxT = np.max(deltaT)
        minT = np.min(deltaT)
        tBins = int( (maxT - minT)/0.0001)
        print('iPDF, tBins, minT, maxT', iPDF, tBins, minT, maxT)
        if tBins == 0:
            continue

        ix = int(iPDF/7)
        iy = iPDF - 7*ix
        axs[ix, iy] = plt.subplot(5, 7, 1+iPDF)
        axs[ix, iy].hist(deltaT, bins=tBins, range=(minT, maxT))
        axs[ix, iy].set_ylabel('at %2.1f eV'%nuMass2)

    
    picName = './spectra/fineSpec/dT2D_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_dataMH%d_fitMH%d_Tmax0.02'%(nuMass1, dist, Ethr, group, dataMH, fitMH)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')



    #    maxNll = np.max(nllVal)
    #    minNll = np.min(nllVal)
    #    nllBins = int( (maxNll - minNll)/1.0)
    #    print('iPDF, nllBins, minNll, maxNll', iPDF, nllBins, minNll, maxNll)
    #    if nllBins == 0:
    #        continue
    #    ix = int(iPDF/7)
    #    iy = iPDF - 7*ix
    #    axs[ix, iy] = plt.subplot(3, 7, 1+iPDF)
    #    axs[ix, iy].hist(nllVal, bins=nllBins, range=(minNll, maxNll))
    #    axs[ix, iy].set_ylabel('at' +nuMass2+ 'eV')



    #picName = './spectra/fineSpec/nll_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_dataMH%d_fitMH%d_Tmax0.02'%(nuMass1, dist, Ethr, group, dataMH, fitMH)
    #plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')



    # plot the nll profile
    crossings = np.zeros(num_datasets) - 1.
    crossings_mask = np.zeros(num_datasets, dtype=bool)
    fig2, axs2 = plt.subplots(figsize=(6, 4))
    for iData in range(num_datasets):
        locMin = np.min(lambdas[iData])
        nllShifted = lambdas[iData] - locMin
        plt.plot(xvals,nllShifted,'-',label='ToyData{}'.format(iData),linewidth=2)
        minNuMass = xvals[np.argmin(lambdas[iData])]
        if minNuMass > 0.7 :
            minDeltaT = np.argmin(lambdas[iData])
            print("EventId %d -> nuMass: %.1f has minimum NLL %.3f"%(iData, minNuMass, locMin/2.))
        
        #if iData%20==0:
        #    print(iData)
        # only fit around the threhsold...
        mask = (nllShifted>=0.)&(nllShifted<7.)
        #print(mask)

        if len(xvals[mask]) > 0:
            p = np.polyfit(xvals[mask],nllShifted[mask],2.)
            crossings[iData] = (-p[1] + np.sqrt( p[1]**2 - 4*(p[0])*(p[2]-2.706) ))/(2*p[0])
            if math.isnan(crossings[iData])==True:
                crossings_mask[iData] = True
            #elif crossings[iData]>1.5:
            #    crossings_mask[iData] = False


    meanSens, medianSens = 0., 0.
    maskedCrossings = np.zeros(num_datasets) - 1.
    maskedCrossings = ma.masked_array(crossings, crossings_mask)
    print(crossings)
    print(crossings_mask)
    print(maskedCrossings)
    meanSens = maskedCrossings.mean()
    medianSens = ma.median(maskedCrossings)

    print(meanSens,'---mean----')
    print(medianSens,'---median----')



    plt.plot(xvals,np.ones(len(xvals))*2.706,'--k')
    plt.text(0.1,2.9,'90% confidence threshold (Wilks\')',fontsize=14)
    #plt.plot(xvals,np.ones(len(xvals))*0.96,'--k')
    #plt.text(0.5,1.2,'90% confidence threshold (toyMC)',fontsize=14)
    plt.axis([0.,3.,0.,30.])
    #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
    plt.ylabel('Test statistic')
    plt.xlabel(r'$m_{\nu}$ (eV)')
    #plt.xlim(0., 0.5)
    #plt.ylim(0, 50)
    picName = './spectra/fineSpec/numass2D_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_dataMH%d_fitMH%d_Tmax0.02'%\
              (nuMass1, dist, Ethr, group, dataMH, fitMH)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    #plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')


    fig3, axs3 = plt.subplots(figsize=(6, 4))
    #plt.hist(crossings, bins=20, range=(0,2))
    nCounts, bins, patches = plt.hist(maskedCrossings,bins=np.linspace(0,3,61))
    plt.xlabel(r"$m_\nu / eV$")
    plt.plot([meanSens, meanSens], [0., np.max(nCounts)],'--k')
    plt.plot([medianSens, medianSens], [0., np.max(nCounts)],':k')
    plt.text(1.5, 0.9*np.max(nCounts), 'mean: %3.2f eV'%meanSens)
    plt.text(1.5, 0.8*np.max(nCounts), 'median: %3.2f eV'%medianSens)


    plt.show()
