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

    num_datasets = 10 #500
    numass_hypotheses = 30
    nStatsPerGroup = 100 #500
    xvals = np.linspace(0, 3.0, numass_hypotheses)
    lambdas = np.zeros((num_datasets, numass_hypotheses))
    bestfitT = np.zeros((num_datasets, numass_hypotheses))
    bestfitN = np.zeros((num_datasets, numass_hypotheses))

    
    for iPDF in range(numass_hypotheses):
        
        groupNum, deltaT, nllVal, nsig = [], [], [], []
        nuMass2 = iPDF * 0.1 #eV
        filename = ns.fitResFileName(modelNum, chaname, dataMH, nuMass1, fitMH, nuMass2, dist, Ethr, group)[1]
        #print(filename)
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
            bestfitT[iData, iPDF] = deltaT[iData]
            bestfitN[iData, iPDF] = nsig[iData]

    

    plt.figure(0)
    for iData in range(num_datasets) :
        plt.plot(xvals, bestfitT[iData], "-", lw=2)

    plt.xlabel(r"$m_\nu / eV$")
    plt.ylabel(r"$\delta T /s$")


    plt.figure(1)
    for iData in range(num_datasets) :
        plt.plot(xvals, bestfitN[iData], "-", lw=2)

    plt.xlabel(r"$m_\nu / eV$")
    plt.ylabel("Nsig")


    plt.figure(2)
    for iData in range(num_datasets) :
        plt.plot(xvals, lambdas[iData] - np.min(lambdas[iData]), "-", lw=2)

    plt.xlabel(r"$m_\nu / eV$")
    plt.ylabel("n.l.l")




    plt.show()
