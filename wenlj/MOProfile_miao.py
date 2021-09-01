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

    dist = float( sys.argv[6] )# unit: kpc
    print('dist: ', dist, ' kpc')

    Ethr = float( sys.argv[7] )# unit: kpc
    print('fit Ethr: ', Ethr, ' MeV')

    group = int( sys.argv[8] )# group number
    print('group number: ', group)




    ## -----> Load Fit Summary Results <----- ##
    #prefix  = "./dataset/%2.1fkpc/TEvisDATA_" %(dist)
    #prefix += "mod%d_cha%d%s_mh%d" %(modelNum, chaname, chaName[chaname], dataMH) 
    #resFilename = prefix + "_data%2.1feV_mh%d_pdf" %(nuMass1, fitMH)
    #suffix = "eV_%2.1fkpc_Ethr%2.1fMeV_group%d_1DFitSummary_Tmax0.02.txt" %(dist, Ethr, group)
    #print(resFilename+suffix)

    num_datasets = 100 #500
    numass_hypotheses = 20
    nStatsPerGroup = 100 #500

    xvalsi_NO = np.linspace(0, 2.0, numass_hypotheses)
    lambdas_NO = np.zeros((num_datasets, numass_hypotheses))
    bestfitT_NO = np.zeros((num_datasets, numass_hypotheses))
    xvalsi_IO = np.linspace(0, 2.0, numass_hypotheses)
    lambdas_IO = np.zeros((num_datasets, numass_hypotheses))
    bestfitT_IO = np.zeros((num_datasets, numass_hypotheses))

    
    for iPDF in range(numass_hypotheses):

        for fitMH in range(1, 3, 1) :
        
            groupNum, deltaT, nllVal, nsig = [], [], [], []
            nuMass2 = iPDF * 0.1 #eV
            filename = ns.fitResFileName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)[1]
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
                if fitMH == 1:
                    lambdas_NO[iData, iPDF] = nllVal[iData] 
                    bestfitT_NO[iData, iPDF] = deltaT[iData] 
                if fitMH == 2:
                    lambdas_IO[iData, iPDF] = nllVal[iData] 
                    bestfitT_IO[iData, iPDF] = deltaT[iData] 


    # specify event id : 
    nuMassArr = np.arange(0, 2.0, 0.1)
    evtId = 1
    plt.plot(nuMassArr, lambdas_NO[evtId] - np.min(lambdas_IO[evtId]), "o-", color="blue", label="NO")
    plt.plot(nuMassArr, lambdas_IO[evtId] - np.min(lambdas_IO[evtId]), "o-", color="orange", label="IO")

    plt.text(0.2, np.mean(lambdas_NO[evtId]), "True MO: Normal", fontsize=14)

    plt.xlabel(r"$m_\nu / eV$")
    plt.ylabel(r"$\Delta \chi^2$")
    plt.legend()
    plt.show()


