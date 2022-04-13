import os, sys
import numpy as np
import matplotlib.pyplot as plt
import math
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

    nuMass2 = float( sys.argv[4] )# unit: eV
    print('nuMass for producing dataset : ', nuMass2, ' eV')
    nuMassName1 = int(nuMass2 * 10)

    MHName = ['NoOsc', 'NH', 'IH']
    dataMH = int( sys.argv[5] )
    print('NMO for dataset : ', dataMH, MHName[dataMH])

    fitMH = int( sys.argv[6] )
    print('NMO for fitting : ', fitMH, MHName[fitMH])

    dist = float( sys.argv[7] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[8] )# unit: MeV
    print('fit Ethr: ', Ethr, ' MeV')

    fitEmax = float( sys.argv[9])   # unit : MeV
    print('fit Emax: ', fitEmax, 'MeV')

    group = int( sys.argv[10] )# group number
    print('group number: ', group)




    ## -----> Load Fit Summary Results <----- ##
    #prefix  = "./dataset/%2.1fkpc/TEvisDATA_" %(dist)
    #prefix += "mod%d_cha%d%s_mh%d" %(modelNum, chaname, chaName[chaname], dataMH) 
    #resFilename = prefix + "_data%2.1feV_mh%d_pdf" %(nuMass1, fitMH)
    #suffix = "eV_%2.1fkpc_Ethr%2.1fMeV_group%d_1DFitSummary_Tmax0.02.txt" %(dist, Ethr, group)
    #print(resFilename+suffix)

    num_datasets = 500
    numass_hypotheses = 20
    nStatsPerGroup = 500
    xvals = np.linspace(0, 2.0, numass_hypotheses)
    #xvals = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9]
    lambdas = np.zeros((num_datasets, numass_hypotheses))
    bestfitT = np.zeros((num_datasets, numass_hypotheses))

    
    resPath = "/junofs/users/miaoyu/supernova/analysis/fitResults/"

    ######## channel eES ##########
    #for iPDF in range(numass_hypotheses):
    #    groupNum, deltaT, nllVal, nsig = [], [], [], []
    #    nuMass1 = iPDF * 0.1 #eV
    #    resFiles = [resPath+"tres_mod%d_cha%d_mh%d%d_nuMass%.1feV.txt"%(modelNum, chaname, dataMH, fitMH, nuMass1), resPath+"res_mod%d_cha%d_mh%d%d_nuMass%.1feV_summary.txt"%(modelNum, chaname, dataMH, fitMH, nuMass1)]
    #    filename = resFiles[1]
    #    print(filename)
    #    if not os.path.exists(filename) :
    #        print("%s does not exist!")
    #        continue
    #    t1, t2, t3, t4 = getNLL(filename, nStatsPerGroup)
    #    groupNum = groupNum + t1
    #    deltaT = deltaT + t2
    #    nllVal = nllVal + t3
    #    nsig   = nsig   + t4

    #    for iData in range(num_datasets):
    #        lambdas[iData, iPDF] = nllVal[iData] * 2


    chaname = 0
    ######### channel pES ##########
    for iPDF in range(numass_hypotheses):
        groupNum, deltaT, nllVal, nsig = [], [], [], []
        nuMass1 = iPDF * 0.1 #eV
        resFiles = [resPath+"res_mod%d_cha%d_mh11_nuMass%.1feV.txt"%(modelNum, chaname, nuMass1), resPath+"res_mod%d_cha%d_mh11_nuMass%.1feV_summary.txt"%(modelNum, chaname, nuMass1)]
        filename = resFiles[1]
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
            lambdas[iData, iPDF] += nllVal[iData] * 2


    # plot the nll profile
    crossings = np.zeros(num_datasets) - 1.
    crossings_mask = np.zeros(num_datasets, dtype=bool)
    fig2, (axs2, axs3) = plt.subplots(1, 2, figsize=(9, 4))
    for iData in range(num_datasets):
        locMin = np.min(lambdas[iData])
        nllShifted = lambdas[iData] - locMin
        axs2.plot(xvals,nllShifted,'-',label='ToyData{}'.format(iData),linewidth=1)
        minNuMass = xvals[np.argmin(lambdas[iData])]
        #if minNuMass > 0.7 :
        #    minDeltaT = np.argmin(lambdas[iData])
        #    print("EventId %d -> nuMass: %.1f has minimum NLL %.3f"%(iData, minNuMass, locMin/2.))
        
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
    #print(crossings)
    #print(crossings_mask)
    #print(maskedCrossings)
    meanSens = maskedCrossings.mean()
    medianSens = ma.median(maskedCrossings)

    print(meanSens,'---mean----')
    print(medianSens,'---median----')



    axs2.plot(xvals,np.ones(len(xvals))*2.706,'--k')
    #plt.text(0.1,2.9,'90% confidence threshold (Wilks\')',fontsize=14)
    #plt.plot(xvals,np.ones(len(xvals))*0.96,'--k')
    #plt.text(0.5,1.2,'90% confidence threshold (toyMC)',fontsize=14)
    #plt.axis([0.,2.,0.,20.])
    #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
    axs2.set_ylabel(r'$\Delta \chi^2$')
    axs2.set_xlabel(r'$m_{\nu}$ [eV]')
    #axs2.set_ylim(0, 10)
    #plt.xlim(0., 0.5)
    #plt.ylim(0, 50)

    #picName = './spectra/fineSpec/numass2D_data%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_dataMH%d_fitMH%d_Tmax0.02'%\
    #          (nuMass1, dist, Ethr, group, dataMH, fitMH)
    #picName = ns.nllProfile2NuMassName(modelNum, chaname, nuType, dataMH, nuMass1, fitMH,  Ethr, group, dist)
    #picName = "nllProfile_stat10000_dataMH%d_fitMH%d_Emax10MeV.pdf"%(dataMH, fitMH)
    #plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    #plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')


    #fig3, axs3 = plt.subplots(figsize=(6, 4))
    #plt.hist(crossings, bins=20, range=(0,2))
    nCounts, bins, patches = plt.hist(maskedCrossings,bins=30, range=(0, 3))
    axs3.set_xlabel(r"$m_\nu$ [eV]")
    axs3.plot([meanSens, meanSens], [0., np.max(nCounts)],'--k')
    axs3.plot([medianSens, medianSens], [0., np.max(nCounts)],':k')
    axs3.text(1.5, 0.9*np.max(nCounts), 'mean: %3.2f eV'%meanSens)
    axs3.text(1.5, 0.8*np.max(nCounts), 'median: %3.2f eV'%medianSens)
    print("statistics error : %.3f" %np.std(maskedCrossings))

    #picName = ns.nuMassUpperLimitName(modelNum, chaname, nuType, dataMH, nuMass1, fitMH, Ethr, group, dist)
    picName = "nuMassUpperLimit_pES_dataMH%d_fitMH%d.pdf"%(dataMH, fitMH)
    plt.tight_layout()
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')

    plt.show()
