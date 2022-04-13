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
    

    fig = plt.figure(figsize=(20, 12))
    spec = gridspec.GridSpec(ncols=8, nrows=4)
    
    ax =  []
    ax0 = fig.add_subplot(spec[0]) ; ax.append(ax0)
    ax1 = fig.add_subplot(spec[1]) ; ax.append(ax1)
    ax2 = fig.add_subplot(spec[2]) ; ax.append(ax2)
    ax3 = fig.add_subplot(spec[3]) ; ax.append(ax3)
    ax4 = fig.add_subplot(spec[4]) ; ax.append(ax4)
    ax5 = fig.add_subplot(spec[5]) ; ax.append(ax5)
    ax6 = fig.add_subplot(spec[6]) ; ax.append(ax6)
    ax7 = fig.add_subplot(spec[7]) ; ax.append(ax7)
    ax8 = fig.add_subplot(spec[8]) ; ax.append(ax8)
    ax9 = fig.add_subplot(spec[9]) ; ax.append(ax9)
    ax10 = fig.add_subplot(spec[10]) ; ax.append(ax10)
    ax11 = fig.add_subplot(spec[11]) ; ax.append(ax11)
    ax12 = fig.add_subplot(spec[12]) ; ax.append(ax12)
    ax13 = fig.add_subplot(spec[13]) ; ax.append(ax13)
    ax14 = fig.add_subplot(spec[14]) ; ax.append(ax14)
    ax15 = fig.add_subplot(spec[15]) ; ax.append(ax15)
    ax16 = fig.add_subplot(spec[16]) ; ax.append(ax16)
    ax17 = fig.add_subplot(spec[17]) ; ax.append(ax17)
    ax18 = fig.add_subplot(spec[18]) ; ax.append(ax18)
    ax19 = fig.add_subplot(spec[19]) ; ax.append(ax19)
    ax20 = fig.add_subplot(spec[20]) ; ax.append(ax20)
    ax21 = fig.add_subplot(spec[21]) ; ax.append(ax21)
    ax22 = fig.add_subplot(spec[22]) ; ax.append(ax22)
    ax23 = fig.add_subplot(spec[23]) ; ax.append(ax23)
    ax24 = fig.add_subplot(spec[24]) ; ax.append(ax24)
    ax25 = fig.add_subplot(spec[25]) ; ax.append(ax25)
    ax26 = fig.add_subplot(spec[26]) ; ax.append(ax26)
    ax27 = fig.add_subplot(spec[27]) ; ax.append(ax27)
    ax28 = fig.add_subplot(spec[28]) ; ax.append(ax28)
    ax29 = fig.add_subplot(spec[29]) ; ax.append(ax29)
    ax30 = fig.add_subplot(spec[30]) ; ax.append(ax30)
    ax31 = fig.add_subplot(spec[31]) ; ax.append(ax31)

    #fig, ax = plt.subplots()
    Nevt = 5000

    test_flag = True
    if test_flag :
        fig1, axx = plt.subplots()

    rightNO, rightIO = [], []

    cc = 0
    tau1_arrNO, tau2_arrNO, A_arrNO, Nsig_arrNO = [], [], [], []
    tau1_arrIO, tau2_arrIO, A_arrIO, Nsig_arrIO = [], [], [], []
    mod_arr = [81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502, 92503,92700,92701,92702,92703,94000,94001,94002,94003]
    for imod in mod_arr:
        train_x, train_y = [[], []], []
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
                    ax[cc].scatter(p0/p1, p2, s=5, alpha=0.3, c="royalblue", label="MO = %d"%imh)
                if imh == 1:
                    ax[cc].scatter(p0/p1, p2, s=5, alpha=0.3, c="crimson",   label="MO = %d"%imh)


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
            
            test_x = np.array(test_x)
            test_x = test_x.T

            pred_y = clf.predict(test_x)
            N = 0
            for n in range(500):
                if pred_y[n] == 1:
                    N += 1

            rightNO.append(N)

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
            
            test_x = np.array(test_x)
            test_x = test_x.T

            pred_y = clf.predict(test_x)
            N = 0
            for n in range(500):
                if pred_y[n] == 2:
                    N += 1

            rightIO.append(N)


        
        #clf = svm.LinearSVC(max_iter=5000)
        #clf.fit(train_x, train_y)

        #clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1)
        #clf.fit(train_x, train_y)

        #print("Iteration times %d"%clf.n_iter_)

            #if imh == 1:
            #    for (i, j, k, l, m) in zip(p0, p1, p2, p3, p4):
            #        if m == 0:
            #            tau1_arrNO.append(i)
            #            tau2_arrNO.append(j)
            #            A_arrNO.append(k)
            #            Nsig_arrNO.append(l)
            #if imh == 2:
            #    for (i, j, k, l, m) in zip(p0, p1, p2, p3, p4):
            #        if m == 0:
            #            tau1_arrIO.append(i)
            #            tau2_arrIO.append(j)
            #            A_arrIO.append(k)
            #            Nsig_arrIO.append(l)


        coef = clf.coef_
        intercept = clf.intercept_
        print("FIT", coef, intercept)
        dx = np.arange(0, 40, 1)
        dy = line(dx, coef, intercept)
        ax[cc].plot(dx, dy, lw=2, color="black")


        ax[cc].set_xlim(0, 30)
        ax[cc].set_ylim(0, 60)
        ax[cc].set_xlabel(r"$\tau_1/\tau_2$", fontsize=9)
        ax[cc].set_ylabel("A", fontsize=9)
        ax[cc].set_title("model %d"%imod, fontsize=9)
        cc += 1


    #tau1_arrNO = np.array(tau1_arrNO)
    #tau2_arrNO = np.array(tau2_arrNO)
    #A_arrNO = np.array(A_arrNO)
    #Nsig_arrNO = np.array(Nsig_arrNO)
    #tau1_arrIO = np.array(tau1_arrIO)
    #tau2_arrIO = np.array(tau2_arrIO)
    #A_arrIO = np.array(A_arrIO)
    #Nsig_arrIO = np.array(Nsig_arrIO)


            #ax[cc].hist(p2, bins=50, range=(0, 0.1),histtype="step", lw=1.5, label="MO = %d"%imh)
        #ax[cc].legend()
        #ax[cc].set_title("Model %d"%imod)
        #ax[cc].set_xlim(0, 40)
    

    #ax.hist(tau1_arrNO, bins=50, histtype="step", color="royalblue")
    #ax.hist(tau1_arrIO, bins=50, histtype="step", color="hotpink")
    #ax.hist(tau2_arrNO, bins=50, histtype="step", color="royalblue")
    #ax.hist(tau2_arrIO, bins=50, histtype="step", color="hotpink")
    #ax.hist(A_arrNO, bins=50, range=(0, 50), histtype="step", color="royalblue")
    #ax.hist(A_arrIO, bins=50, range=(0, 50), histtype="step", color="hotpink")
    #ax.hist(Nsig_arrNO, bins=50, histtype="step", color="royalblue")
    #ax.hist(Nsig_arrIO, bins=50, histtype="step", color="hotpink")
    #ax.hist(tau1_arrNO/tau2_arrNO, bins=50, range=(0, 10), histtype="step", color="royalblue")
    #ax.hist(tau1_arrIO/tau2_arrIO, bins=50, range=(0, 10), histtype="step", color="hotpink")

    #ax.scatter(tau1_arrNO/tau2_arrNO, A_arrNO, s=1, c="royalblue", label="NO")
    #ax.scatter(tau1_arrIO/tau2_arrIO, A_arrIO, s=1, c="hotpink", alpha=0.1, label="IO")
    ##ax.scatter(tau2_arrNO, tau1_arrNO,  s=1, c="royalblue")
    ##ax.scatter(tau2_arrIO, tau1_arrIO,  s=1, c="crimson")

    #ax.legend(prop={'size':14})
    #ax.set_xlabel(r"$\tau_1/\tau_2$", fontsize=14)
    #ax.set_ylabel("A", fontsize=14)

    #ax.set_xlim(0, 20)
    #ax.set_ylim(0, 60)

    rightNO = np.array(rightNO)
    rightIO = np.array(rightIO)

    mod_arr = np.arange(0, 32, 1)
    axx.plot(mod_arr, rightNO/500., "^-", color="crimson", label="Normal Ordering")
    axx.plot(mod_arr, rightIO/500., "v-", color="royalblue",   label="Inverted Ordering")
    axx.legend(prop={"size":13})
    axx.set_xlabel("Model No.", fontsize=13)
    axx.set_ylabel("efficiency")


    plt.tight_layout()
    #plt.savefig("scatter_2kpc.pdf")
    plt.show()








if __name__ == "__main__" :
    draw_scatter()
