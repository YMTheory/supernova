import uproot as up
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys, os

from sklearn.linear_model import SGDClassifier
from sklearn import svm
from sklearn.neural_network import MLPClassifier

import ROOT

def line(x, coef, intercept):
    return (-(x * coef[0, 0]) - intercept[0]) / coef[0, 1]

def draw_scatter():
    
    fig, ax = plt.subplots()

    Nevt = 5000

    test_flag = True
    if test_flag :
        fig1, axx = plt.subplots()

    cc = 0
    tau1_arrNO, tau2_arrNO, A_arrNO, Nsig_arrNO = [], [], [], []
    tau1_arrIO, tau2_arrIO, A_arrIO, Nsig_arrIO = [], [], [], []
    mod_arr = [81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502, 92503,92700,92701,92702,92703,94000,94001,94002,94003]
    train_x, train_y = [[], []], []
    for imod in mod_arr:
        nData  = 0
        for imh in range(1, 3):
        #for imod in [82503]:
            for no in range(10, 11, 1):
                filename = "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod%d_MH%d_train%d.root"%(imod, imh, no)
                if filename == "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod82703_MH1_train0.root":
                    continue
                if filename == "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod82703_MH2_train3.root":
                    continue
                if filename == "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod91122_MH2_train1.root":
                    continue
                print(filename)
                ff = up.open(filename)
                p0 = ff["fit"]["tau1"].array()
                p1 = ff["fit"]["tau2"].array()
                p2 = ff["fit"]["A"].array()
                p3 = ff["fit"]["Nsig"].array()
                p4 = ff['fit']['status'].array()

                if imh == 2:
                    ax.scatter(p0/p1, p2, s=5, alpha=0.3, c="royalblue", label="MO = %d"%imh)
                if imh == 1:
                    ax.scatter(p0/p1, p2, s=5, alpha=0.3, c="crimson",   label="MO = %d"%imh)


                if imh == 1:
                    for n in range(Nevt):
                        if 0 < p0[n] / p1[n] < 60 and 0 < p2[n] < 30 :
                            train_x[0].append(p0[n]/p1[n])
                            train_x[1].append(p2[n])
                            train_y.append(1)
                        #train_x[n, 0] = p0[n] / p1[n]
                        #train_x[n, 1] = p2[n]
                        #train_y[n] = 1
                if imh == 2:
                    for n in range(Nevt):
                        if 0 < p0[n] / p1[n] < 60 and 0 < p2[n] < 30 :
                            train_x[0].append(p0[n]/p1[n])
                            train_x[1].append(p2[n])
                            train_y.append(2)
                        #train_x[n+Nevt, 0] = p0[n] / p1[n]
                        #train_x[n+Nevt, 1] = p2[n]
                        #train_y[n+Nevt] = 2


    train_x = np.array(train_x)
    train_x = train_x.T
    train_y = np.array(train_y)
        
    clf = SGDClassifier(loss="hinge", penalty="l2", max_iter=5000)
    clf.fit(train_x, train_y)

    if test_flag :
        test_x, test_y = [[], []], []

        for imod in mod_arr:

            ###### Normal Mass Ordering ####
            filename = "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod%d_MH%d_test10.root"%(imod, 1)
            print(filename)
            ff = up.open(filename)
            p0 = ff["fit"]["tau1"].array()
            p1 = ff["fit"]["tau2"].array()
            p2 = ff["fit"]["A"].array()
            p3 = ff["fit"]["Nsig"].array()
            p4 = ff['fit']['status'].array()

            for n in range(500):
                test_x[0].append(p0[n]/p1[n])
                test_x[1].append(p2[n])
                test_y.append(1)
            
        test_x = np.array(test_x)
        test_x = test_x.T

        pred_y = clf.predict(test_x)
        N = 0
        for n in range(len(test_y)):
            if pred_y[n] == 1:
                N += 1
    
        print("NO efficiency = %.3f"%(N/float(len(test_y))) )

        test_x, test_y = [[], []], []
        for imod in mod_arr:
            ###### Inverted Mass Ordering ####
            test_x, test_y = [[], []], []
            filename = "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod%d_MH%d_test10.root"%(imod, 2)
            print(filename)
            ff = up.open(filename)
            p0 = ff["fit"]["tau1"].array()
            p1 = ff["fit"]["tau2"].array()
            p2 = ff["fit"]["A"].array()
            p3 = ff["fit"]["Nsig"].array()
            p4 = ff['fit']['status'].array()

            for n in range(500):
                test_x[0].append(p0[n]/p1[n])
                test_x[1].append(p2[n])
                test_y.append(2)
            
        test_x = np.array(test_x)
        test_x = test_x.T

        pred_y = clf.predict(test_x)
        N = 0
        for n in range(len(test_y)):
            if pred_y[n] == 2:
                N += 1


        print("IO efficiency = %.3f"%(N/float(len(test_y))) )

        

    coef = clf.coef_
    intercept = clf.intercept_
    print("FIT", coef, intercept)
    dx = np.arange(0, 40, 1)
    dy = line(dx, coef, intercept)
    ax.plot(dx, dy, lw=2, color="black")


    ax.set_xlim(0, 30)
    ax.set_ylim(0, 60)
    ax.set_xlabel(r"$\tau_1/\tau_2$", fontsize=9)
    ax.set_ylabel("A", fontsize=9)
    cc += 1



    plt.tight_layout()
    #plt.savefig("scatter_2kpc.pdf")
    plt.show()








if __name__ == "__main__" :
    draw_scatter()
