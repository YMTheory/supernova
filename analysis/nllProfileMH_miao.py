import numpy as np
import matplotlib.pyplot as plt
import sys, os
import math
import numpy.ma as ma



def getValue(strline):
    values = strline.split()
    x1 = int(values[0])
    x2 = float(values[1])
    x3 = float(values[2])
    x4 = float(values[3])
    x5 = int(values[4])
    x6 = int(values[5])
    return (x1, x2, x3, x4, x5, x6)

def func(x, a0, t0, a1):
    return a0*np.power(x-t0, 2) + a1

def getNLL(filename, num_datasets):
    groupNum, deltaT, nllVal, nsig = [], [], [], []
    trueStat, fitStatus = [], []
    
    if os.path.exists(filename) == False:
        return groupNum, deltaT, nllVal, nsig, trueStat, fitStatus

    file = open(filename, 'r')
    lines = file.readlines()
    file.close()

    if num_datasets > len(lines):
        print(filename)
        print('Error! num_datasets is larger than the number of rows in txt file!')
        return groupNum, deltaT, nllVal, nsig, trueStat, fitStatus

    for k in range(num_datasets):
        (num, dT, val, nEvt, nMC, fitS) = getValue( lines[k] )
        groupNum.append(num)
        deltaT.append(dT)
        nllVal.append(val)
        nsig.append(nEvt)
        trueStat.append(nMC)
        fitStatus.append(fitS)
        if fitS!=0:
            print('fitS!=0: num, k, status', num, k, fitS)
    return groupNum, deltaT, nllVal, nsig, trueStat, fitStatus




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

    dist = float( sys.argv[6] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[7] )# unit: kpc
    print('fit Ethr: ', Ethr, ' MeV')

    Emax = float( sys.argv[8] )# unit: kpc
    print('fit Emax: ', Emax, ' MeV')

    group = int( sys.argv[9] )# group number
    print('group number: ', group)

    
    Tmax = 0.02
    

    num_datasets = 500
    numass_hypotheses = 20

    xvals_NO = np.linspace(0, 2.0, numass_hypotheses)
    xvals_IO = np.linspace(0, 2.0, numass_hypotheses)

    lambdas_NO = np.zeros((num_datasets, numass_hypotheses)) 
    lambdas_IO = np.zeros((num_datasets, numass_hypotheses))

    bestfitT_NO = np.zeros((num_datasets, numass_hypotheses))
    bestfitT_IO = np.zeros((num_datasets, numass_hypotheses))

    nllStatus_NO = np.zeros((num_datasets, numass_hypotheses))
    nllStatus_IO = np.zeros((num_datasets, numass_hypotheses))

    deltaChi2 = np.zeros((num_datasets, numass_hypotheses))

    fig, axs = plt.subplots(3, 7, figsize=(30, 5))

    # Loading fitting results from output txt
    resPath = "/junofs/users/miaoyu/supernova/analysis/fitResults/"
    for iPDF in range(numass_hypotheses):

        groupNum_NO, deltaT_NO, nllVal_NO, nsig_NO = [], [], [], []
        groupNum_IO, deltaT_IO, nllVal_IO, nsig_IO = [], [], [], []
        nStat_NO, status_NO = [], []
        nStat_IO, status_IO = [], []
        
        for imh in range(1, 3, 1):

            nuMass2 = iPDF * 0.1

            #filename = ns.fitResFileName(modelNum, chaname, dataMH, nuMass1, imh, nuMass2, dist, Ethr, Emax, group)[1]
            resFiles = [resPath+"testres_mod%d_cha%d_mh%d%d_nuMass%.1feV.txt"%(modelNum, chaname, dataMH, imh, nuMass1), resPath+"testres_mod%d_cha%d_mh%d%d_nuMass%.1feV_summary.txt"%(modelNum, chaname, dataMH, imh, nuMass1)]
            filename = resFiles[1]
            print(filename)
            t1, t2, t3, t4, t5, t6 = getNLL(filename, num_datasets)

            if imh == 1:   #NO
                groupNum_NO = groupNum_NO + t1
                deltaT_NO   = deltaT_NO + t2
                nllVal_NO   = nllVal_NO + t3
                nsig_NO     = nsig_NO   + t4
                nStat_NO    = nStat_NO  + t5
                status_NO   = status_NO  + t6

            if imh == 2:   #IO
                groupNum_IO = groupNum_IO + t1
                deltaT_IO   = deltaT_IO + t2
                nllVal_IO   = nllVal_IO + t3
                nsig_IO     = nsig_IO   + t4
                nStat_IO    = nStat_IO  + t5
                status_IO   = status_IO  + t6
        
        for iData in range(num_datasets):
            lambdas_NO[iData, iPDF] = nllVal_NO[iData] * 2 # chi2 ~ 2* (NLL(x) - NLL_loc_min)
            lambdas_IO[iData, iPDF] = nllVal_IO[iData] * 2 # chi2 ~ 2* (NLL(x) - NLL_loc_min)
            nllStatus_NO[iData, iPDF] = status_NO[iData]
            nllStatus_IO[iData, iPDF] = status_IO[iData]
            bestfitT_NO[iData, iPDF] = deltaT_NO[iData]
            bestfitT_IO[iData, iPDF] = deltaT_IO[iData]

            if dataMH == 1:
                deltaChi2[iData, iPDF] = 2*(nllVal_IO[iData] - nllVal_NO[iData])

            if dataMH==2:
                deltaChi2[iData, iPDF] = 2*(nllVal_NO[iData] - nllVal_IO[iData])

            
            #print(iData, iPDF, lambdas_NO[iData, iPDF]/2., lambdas_IO[iData, iPDF]/2.)

        maxT_NO = np.max(deltaT_NO)
        maxT_IO = np.max(deltaT_IO)

        minT_NO = np.min(deltaT_NO)
        minT_IO = np.min(deltaT_IO)

        maxT = maxT_NO
        if maxT<maxT_IO:
            maxT = maxT_IO

        minT = minT_NO
        if minT>minT_IO:
            minT = minT_IO

        maxT = maxT*1.1
        minT = minT*0.9
        tBins = int( (maxT - minT)/0.0001 )

        if tBins == 0:
            continue


        ix = int(iPDF/7)
        iy = iPDF - 7*ix

        axs[ix, iy] = plt.subplot(3, 7, 1+iPDF)
        axs[ix, iy].hist(deltaT_IO, bins=tBins, range=(minT, maxT), color='r',\
              histtype='step')
        axs[ix, iy].hist(deltaT_NO, bins=tBins, range=(minT, maxT), facecolor='b',\
              alpha=0.75)
        axs[ix, iy].set_ylabel('at %2.1f eV'%(nuMass2))

    picName  = './fitResults/dT_data%2.1feV_mh%d'%(nuMass1, dataMH)
    picName += '_%2.1fkpc_Ethr%2.1fMeV_Tmax%3.2f_group%d'%(dist, Ethr, Tmax, group)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')


    ######################################################

    chi2Min, chi2Max = -30., 30.
    if dist == 3.0:
        if dataMH == 1:
            chi2Min, chi2Max = 0., 80.
        if dataMH == 2:
            chi2Min, chi2Max = 0., 250.

    if dist == 5.0:
        if dataMH == 1:
            chi2Min, chi2Max = -25., 50.
        if dataMH == 2:
            chi2Min, chi2Max = -25., 150.

    if dist == 10.0:
        if dataMH == 1:
            chi2Min, chi2Max = -25., 50.
        if dataMH == 2:
            chi2Min, chi2Max = -25., 150.

    ######################################################

    meanDchisq, varDchisq, stdDchisq = [], [], []
    fig, axs = plt.subplots(3, 7, figsize=(30, 5))
    for iPDF in range(numass_hypotheses):

        # remove the bad fits
        #mask_NO = (nllStatus_NO[:,iPDF]!=0)
        #mask_IO = (nllStatus_IO[:,iPDF]!=0)
        #mask    = (nllStatus_IO[:,iPDF]!=0) | (nllStatus_NO[:,iPDF]!=0)
        mask_NO = (nllStatus_NO[:,iPDF]==0)
        mask_IO = (nllStatus_IO[:,iPDF]==0)
        mask    = (nllStatus_IO[:,iPDF]==0) | (nllStatus_NO[:,iPDF]==0)

        maskedDchisq = ma.masked_array(deltaChi2[:,iPDF], mask)
        #if iPDF==0:
        #    print('test mask_NO: ', mask_NO)
        #    print('test mask_IO: ', mask_IO)
        #    print('test mask: ', mask)
        #    print('test deltaChi2 original: ', deltaChi2[:,iPDF])
        #    print('test deltaChi2 masked: ', deltaChi2[mask][iPDF])
        #    print('maskedDchisq: ', maskedDchisq)

        axs[ix, iy] = plt.subplot(3, 7, 1+iPDF)
        axs[ix, iy].hist(maskedDchisq, bins=20, range=(chi2Min, chi2Max), \
              facecolor='g', alpha=0.75)

        #deltachisq = np.array(deltaChi2[:, iPDF])
        meanDchisq.append( maskedDchisq.mean() )
        varDchisq.append(  maskedDchisq.var() )
        stdDchisq.append(  maskedDchisq.std() )

    print('meanDchisq: ', meanDchisq)
    print('varDchisq: ', varDchisq)
    print('stdDchisq: ', stdDchisq)
    picName  = './fitResults/deltaChi2_data%2.1feV_mh%d'%(nuMass1, dataMH)
    picName += '_%2.1fkpc_Ethr%2.1fMeV_Tmax%3.2f_group%d'%(dist, Ethr, Tmax, group)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')

    #deltachisq = np.array(deltachisq)
    meanDchisq = np.array(meanDchisq)
    stdDchisq  = np.array(stdDchisq)

    fig, ax = plt.subplots(1,1)
    if dataMH == 1:
        labelName = 'True: NO'

    if dataMH == 2:
        labelName = 'True: IO'
    ax.plot( xvals_NO, meanDchisq, 'bs-', label=labelName )
    ax.fill_between( xvals_NO, meanDchisq + stdDchisq, meanDchisq - stdDchisq, \
                     facecolor='blue', alpha=0.5)
    #nCounts, bins, patches = plt.hist(DeltaChisq,bins=np.linspace(0,chi2Max,100))
    ax.set_ylabel(r'$\Delta\chi^2$')
    ax.set_xlabel(r'$m_{\nu}$ (eV)')
    picName = './fitResults/deltaChisqVSmNu_data%2.1f_mh%d_%2.1fkpc_group%d'\
              %(nuMass1, dataMH, dist, group)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    #plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')

    #fig3, ax3 = plt.subplots(1,1)

    ######################################################
    # plot the nll profile
    crossings_NO = np.zeros(num_datasets) - 1.
    crossings_IO = np.zeros(num_datasets) - 1.

    DeltaChisq = np.zeros(num_datasets) - 1.

    crossings_mask_NO = np.zeros(num_datasets,dtype=bool)
    crossings_mask_IO = np.zeros(num_datasets,dtype=bool)

    ################### analysis
    locMin = 999.
    nllShifted_NO = np.zeros((num_datasets, numass_hypotheses))
    nllShifted_IO = np.zeros((num_datasets, numass_hypotheses))

    for iData in range(num_datasets):

        locMin_NO = np.min(lambdas_NO[iData])
        locMin_IO = np.min(lambdas_IO[iData])
        #print(locMin_NO/2., locMin_IO/2.)
        if dataMH == 1:
            locMin = locMin_NO
            DeltaChisq[iData] = locMin_IO - locMin_NO

        if dataMH == 2:
            locMin = locMin_IO
            DeltaChisq[iData] = locMin_NO - locMin_IO

        nllShifted_NO[iData] = lambdas_NO[iData] - locMin
        nllShifted_IO[iData] = lambdas_IO[iData] - locMin
        #print('locMin_NO, locMin_NO: %10.3f  %10.3f'%(locMin_NO, locMin_IO) )

        if iData%20==0:
            print(iData)

        # ##################################
        # check if the slope switches from positive to negative
        #isBadCurve = False
        #if nllShifted_NO[iData][0] >2.:
        #    print('nllShifted_NO[iData][0] >2.: %d'%(iData))
        #    print(nllShifted_NO[iData])
        #    print(lambdas_NO[iData])
        #    print(lambdas_IO[iData]-lambdas_NO[iData])
            
        #for kk in range(2, numass_hypotheses-1):
        #    if dataMH == 1:
        #        slopePre = nllShifted_NO[iData][kk]-nllShifted_NO[iData][kk-1]
        #        slopePost = nllShifted_NO[iData][kk]-nllShifted_NO[iData][kk+1]
        #        if slopePre>0. and slopePost>0.:
        #            isBadCurve = True
        #            #print('nll profile not good, iData, kk: %d, %d, %10.3f, %10.3f'\
        #            #      %(iData, kk, slopePre, slopePost))
        #    
        #    if dataMH == 2:
        #        slopePre = nllShifted_IO[iData][kk]-nllShifted_IO[iData][kk-1]
        #        slopePost = nllShifted_IO[iData][kk]-nllShifted_IO[iData][kk+1]
        #        if slopePre>0. and slopePost>0.:
        #            isBadCurve = True
        #            #print('nll profile not good, iData, kk: %d, %d, %10.3f, %10.3f'\
        #            #      %(iData, kk, slopePre, slopePost))

        #if isBadCurve==True:
        #    print('nll profile not good, iData: %d', iData)
        #    #continue

        ########################
        # the normal ordering case
        # only fit around the threhsold...
        if dataMH == 1:
            mask = (nllStatus_NO[iData,:]!=0)
            maskedNLL = ma.masked_array(nllShifted_NO[iData,:], mask)
            maskedX = ma.masked_array(xvals_NO, mask)

            if nllShifted_NO[iData,0]>10.:
                mask = np.ones(numass_hypotheses, dtype=bool)
                print('mask iData: ', iData, mask)
            else:
                mask = (nllShifted_NO[iData]<0.) | (nllShifted_NO[iData]>7.)

            if np.count_nonzero(mask)!=20:
                p = np.polyfit(maskedX, maskedNLL, 2.)
                #crossings_NO[iData] = (-p[1] + np.sqrt( p[1]**2 - 4*(p[0])*(p[2]-0.88) ))/(2*p[0])
                crossings_NO[iData] = (-p[1] + np.sqrt( p[1]**2 - 4*(p[0])*(p[2]-2.706) ))/(2*p[0])
                if math.isnan(crossings_NO[iData])==True:
                    print('2 crossings_NO[iData]==math.nan')
                    crossings_mask_NO[iData] = True

        ########################
        # the inverted ordering case
        if dataMH == 2:
            mask = (nllStatus_IO[iData,:]!=0)
            maskedNLL = ma.masked_array(nllShifted_IO[iData,:], mask)
            maskedX = ma.masked_array(xvals_IO, mask)

            if nllShifted_IO[iData,0]>10.:
                mask = np.ones(numass_hypotheses, dtype=bool)
                print('mask iData: ', iData, mask)
            else:
                mask = (nllShifted_IO[iData]<0.) | (nllShifted_IO[iData]>7.)

            if np.count_nonzero(mask)!=20:
                p = np.polyfit(maskedX, maskedNLL, 2.)
                #crossings_IO[iData] = (-p[1] + np.sqrt( p[1]**2 - 4*(p[0])*(p[2]-0.96) ))/(2*p[0])
                crossings_IO[iData] = (-p[1] + np.sqrt( p[1]**2 - 4*(p[0])*(p[2]-2.706) ))/(2*p[0])
                if math.isnan(crossings_IO[iData])==True:
                    print('2 crossings_IO[iData]==math.nan')
                    crossings_mask_IO[iData] = True

    ################### plot
    fig2 = plt.subplots(2, 1, figsize=(3.6, 5.4))
    plt.subplots_adjust(left=0.12, top= 0.97, right = 0.97, bottom = 0.14, wspace = 0.01, hspace = 0.01)

    axs = plt.subplot(2, 1, 1)
    for iData in range(num_datasets):
        if dataMH==2:
            mask = (nllStatus_NO[iData,:]!=0)
            maskedNLL = ma.masked_array(nllShifted_NO[iData,:], mask)
            maskedX = ma.masked_array(xvals_NO, mask)

        if dataMH==1:
            mask = (nllStatus_IO[iData,:]!=0)
            maskedNLL = ma.masked_array(nllShifted_IO[iData,:], mask)
            maskedX = ma.masked_array(xvals_IO, mask)

        if iData==0:
            print('mask: ', mask)
            print('maskedX: ', maskedX)
            print('maskedNLL: ', maskedNLL)
        plt.plot(maskedX,maskedNLL,'-',label='ToyData{}'.format(iData),linewidth=0.5)

    #plt.axis([0.,2.,11.,100.])
    plotYmin, plotYmax = 2., 100.
    if dist==5. or dist==3.:
        plotYmin, plotYmax = 2., 500.

    plt.axis([0.,2.,plotYmin,plotYmax])
    axs.set_yscale('log')
    for label in axs.get_xticklabels():
        label.set_visible(False)

    axs = plt.subplot(2, 1, 2)
    for iData in range(num_datasets):
        if dataMH==2:
            mask = (nllStatus_IO[iData,:]!=0)
            if nllShifted_IO[iData,0]>10.:
                mask = np.ones(numass_hypotheses, dtype=bool)
                print('mask iData: ', iData, mask)

            maskedNLL = ma.masked_array(nllShifted_IO[iData,:], mask)
            maskedX = ma.masked_array(xvals_IO, mask)

        if dataMH==1:
            mask = (nllStatus_NO[iData,:]!=0)
            if nllShifted_NO[iData,0]>10.:
                mask = np.ones(numass_hypotheses, dtype=bool)
                print('mask iData: ', iData, mask)

            maskedNLL = ma.masked_array(nllShifted_NO[iData,:], mask)
            maskedX = ma.masked_array(xvals_NO, mask)

        if np.count_nonzero(mask)!=20:
            plt.plot(maskedX,maskedNLL,'--',label='ToyData{}'.format(iData),linewidth=0.5)

    plt.plot(xvals_NO,np.ones(len(xvals_NO))*2.706,'--k')
    plt.text(0.5,2.9,'90% confidence threshold (Wilks\')',fontsize=10)
    plt.axis([0.,2.,0.,20.])
    plt.ylabel('Test statistic')
    #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
    plt.xlabel(r'$m_{\nu}$ (eV)')
    picName = './fitResults/numass_data%2.1feV_mh%d_%2.1fkpc_Ethr%2.1fMeV_group%d'%\
              (nuMass1, dataMH, dist, Ethr, group)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    #plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')
    #plt.savefig('./spectra/numass_data%2.1feV_%2.1fkpc_group%d.png'%(nuMass1, dist, group), dpi=400, bbox_inches='tight')

    meanSens, medianSens = 0., 0.
    maskedCrossings = np.zeros(num_datasets) - 1.
    if dataMH == 1:
        print('crossings_NO: ', crossings_NO)
        maskedCrossings = ma.masked_array(crossings_NO, crossings_mask_NO)

    if dataMH == 2:
        print('crossings_IO: ', crossings_IO)
        #print('crossings_mask_IO: ', crossings_mask_IO)
        maskedCrossings = ma.masked_array(crossings_IO, crossings_mask_IO)
        #print('maskedCrossings: ', maskedCrossings)

    meanSens = maskedCrossings.mean()
    medianSens = ma.median(maskedCrossings)

    print(meanSens,'---mean----')
    print(medianSens,'---median----')

    fig3, ax3 = plt.subplots(1,1)
    nCounts, bins, patches = plt.hist(maskedCrossings,bins=np.linspace(0,2,41))
    nCounts = np.array(nCounts)
    ax3.set_ylabel('Counts ({} trials total)'.format(num_datasets))
    ax3.set_xlabel('neutrino mass (eV) at 90% confidence limit')
    plt.plot([meanSens, meanSens], [0., np.max(nCounts)],'--k')
    plt.plot([medianSens, medianSens], [0., np.max(nCounts)],':k')
    plt.text(1.5, 0.9*np.max(nCounts), 'mean: %3.2f eV'%meanSens)
    plt.text(1.5, 0.8*np.max(nCounts), 'median: %3.2f eV'%medianSens)
    #plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
    picName = './fitResults/nuMassSens_data%2.1f_mh%d_%2.1fkpc_group%d'\
              %(nuMass1, dataMH, dist, group)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    #plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')

    fig4, ax4 = plt.subplots(1,1)
    print("========>")
    print(chi2Max)
    print(DeltaChisq)
    meanNMOSens = DeltaChisq.mean()
    medianNMOSens = ma.median(DeltaChisq)
    nCounts, bins, patches = plt.hist(DeltaChisq,bins=np.linspace(DeltaChisq.min()-10, DeltaChisq.max()+10,100), color="blue", edgecolor="black")
    nCounts = np.array(nCounts)
    plt.plot([medianNMOSens, medianNMOSens], [0, np.max(nCounts)], "--k")
    plt.text(medianNMOSens+4, 0.9*np.max(nCounts), r"median $\Delta \chi^2$: %3.2f "%medianNMOSens )
    ax4.set_ylabel('Counts ({} trials total)'.format(num_datasets))
    if dataMH == 2:
        ax4.set_xlabel(r'$\chi^2(NO) - \chi^2(IO)$')
    if dataMH == 1:
        ax4.set_xlabel(r'$\chi^2(IO) - \chi^2(NO)$')
    #picName = './spectra/fineSpec/NMOSens_data%2.1f_mh%d_%2.1fkpc_group%d'\
    #          %(nuMass1, dataMH, dist, group)
    picName = "./fitResults/NMOSens_dataMH%d_Emax10MeV.pdf"%(dataMH)
    plt.savefig(picName+'.pdf', dpi=400, bbox_inches='tight')
    #plt.savefig(picName+'.png', dpi=400, bbox_inches='tight')

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











