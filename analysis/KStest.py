import numpy as np
import matplotlib.pyplot as plt
import ROOT
import uproot as up


bin_width = 0.001   # 1 ms
Tmin      = -0.02
Tmax      = 0.1


def cdf(Tstart, Tend, x, y):
    cml = 0
    for i, j in zip(x, y):
        if Tstart < i <= Tend:
            cml += j

    return cml


def deltaCDF(x1, x2):
    deltax = np.abs(x1 - x2)
    return np.max(deltax)



def loadData(filename):
    ff = up.open(filename)
    evtID = ff["evt"]["evtID"].array()
    time = ff["evt"]["time"].array()
    return evtID, time



def projection(ff, hh, cha):
    f1 = ROOT.TFile(ff, "read")
    h2d = f1.Get(hh)
    Emin = 0.1
    if cha == 0:
        Emax = 3
    if cha == 1:
        Emax = 80

    binEmin = h2d.GetYaxis().FindBin(Emin)
    binEmax = h2d.GetYaxis().FindBin(Emax)
    Tmin, Tmax = -0.020, 0.1
    bin1 = h2d.GetXaxis().FindBin(Tmin)
    bin2 = h2d.GetXaxis().FindBin(Tmax)

    binCent, binCont = [], []
    for i in range(bin1, bin2, 1):
        binCent.append(h2d.GetXaxis().GetBinCenter(i))
        tmp_count = 0
        for j in range(binEmin, binEmax, 1):
            tmp_count += h2d.GetBinContent(i, j) 
        binCont.append(tmp_count)

    binCent = np.array(binCent)
    binCont = np.array(binCont)
    return binCent, binCont





def main():

    path = "../production/PDFs/2kpc/evistSpec_mod"
    
    ich = 1
    mod = [82702]
    #mod = [ 81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502,92503,92700,92701,92702,92703,94000,94001,94002,94003]
    lstyles = ["-", "--"]
    colors = ["black", "gray", "blue", "red", "orange", "darkviolet", "green", "salmon", "gold", "cyan", "dodgerblue", "aquamarine", "rosybrown", "burlywood", "royalblue", "mediumpurple", \
              "hotpink", "sienna", "lightgreen", "deepskyblue", "violet", "olive", "indigo", "aquamarine", "peru", "darkslategrey", "springgreen", "maroon", "goldenrod", "teal", "steelblue", "bisque"]

    cfd_arr_NO, cfd_arr_IO = [], []

    fig, ax = plt.subplots()
    cid = 0
    for imod in mod:
        for imh in [1, 2]:

            filename = path + str(imod) + "_cha" + str(ich) + ".root"
            histname = "hET_mod"+str(imod)+"_cha"+str(ich)+"_mh"+str(imh)
            print(filename, histname)
            
            binCent, binCont = projection(filename, histname, ich)

            #ax.plot(binCent*1000, binCont, "-",  lw=2.0, label="mod%d_cha%d_mh%d"%(imod, ich, imh))
            
            if imh == 1:   # NO
                binCent1, binCont1 = [], []
                N0 = cdf(Tmin, Tmax, binCent, binCont)
                for t in np.arange(Tmin, Tmax, bin_width):
                    N = cdf(Tmin, t, binCent, binCont)
                    binCent1.append((t-Tmin) / (Tmax - Tmin))
                    binCont1.append(N / N0)
                binCent1 = np.array(binCent1)
                binCont1 = np.array(binCont1)
                cfd_arr_NO.append(binCont1)
                #ax.plot(binCent1, binCont1/N0, lstyles[imh-1],  lw=2.0, color=colors[cid], label="mod%d_cha%d_mh%d"%(imod, ich, imh))
            if imh == 2:   # IO
                binCent2, binCont2 = [], []
                N0 = cdf(Tmin, Tmax, binCent, binCont)
                for t in np.arange(Tmin, Tmax, bin_width):
                    N = cdf(Tmin, t, binCent, binCont)
                    binCent2.append((t-Tmin) / (Tmax - Tmin))
                    binCont2.append(N / N0)
                binCent2 = np.array(binCent2)
                binCont2 = np.array(binCont2)
                cfd_arr_IO.append(binCont2)
                #ax.plot(binCent1, binCont1/N0, lstyles[imh-1],  lw=2.0, color=colors[cid], label="mod%d_cha%d_mh%d"%(imod, ich, imh))
        cid += 1
        print("****************************************"+str(imod)+"**********************************************")

    #ax.set_xlabel("time [ms]", fontsize=13)
    #ax.set_ylabel("eES events in JUNO per bin", fontsize=13)
    #ax.set_xlabel("x", fontsize=13)
    #ax.set_ylabel("K(x)", fontsize=13)

    #ax.text(0.1, 0.4, "IO", fontsize=13)
    #ax.text(0.6, 0.2, "NO", fontsize=13)

    ##ax.legend(loc="upper left", prop={'size':16}, ncol=2)
    ##ax.legend()

    #plt.tight_layout()
    #plt.savefig("CDF_models.pdf")
    #plt.show()



    ############# calculate distance differences ##########
    #dmin_NO, dmin_IO = [], []

    ### minimum test:
    #for i in range(32):
    #    dmin1NO = 1000000
    #    dmin1IO = 1000000
    #    for j in range(32):
    #        tmpd1NO = deltaCDF(cfd_arr_NO[i], cfd_arr_IO[j])
    #        tmpd1IO = deltaCDF(cfd_arr_IO[i], cfd_arr_NO[j])
    #        if tmpd1NO < dmin1NO :
    #            dmin1NO = tmpd1NO
    #        if tmpd1IO < dmin1IO:
    #            dmin1IO = tmpd1IO


    #    dmin2NO = 1000000
    #    dmin2IO = 1000000
    #    for k in range(32):
    #        if k == i:
    #            continue
    #        tmpd2NO = deltaCDF(cfd_arr_NO[i], cfd_arr_NO[k])
    #        tmpd2IO = deltaCDF(cfd_arr_IO[i], cfd_arr_IO[k])
    #        if tmpd2NO < dmin2NO :
    #            dmin2NO = tmpd2NO
    #        if tmpd2IO < dmin2IO :
    #            dmin2IO = tmpd2IO


    #    dmin_NO.append(dmin1NO - dmin2NO)
    #    dmin_IO.append(dmin1IO - dmin2IO)

    #
    #    print(mod[i], dmin_NO[-1], dmin_IO[-1])


    #dave_NO, dave_IO = [], []

    #### average test:
    #for i in range(32):
    #    dave1NO = 0
    #    dave1IO = 0
    #    for j in range(32):
    #        tmpd1NO = deltaCDF(cfd_arr_NO[i], cfd_arr_IO[j])
    #        tmpd1IO = deltaCDF(cfd_arr_IO[i], cfd_arr_NO[j])
    #        dave1NO += tmpd1NO
    #        dave1IO += tmpd1IO

    #    dave1NO /= 32.
    #    dave1IO /= 32.

    #    dave2NO = 0
    #    dave2IO = 0
    #    for k in range(32):
    #        if k == i:
    #            continue
    #        tmpd2NO = deltaCDF(cfd_arr_NO[i], cfd_arr_NO[k])
    #        tmpd2IO = deltaCDF(cfd_arr_IO[i], cfd_arr_IO[k])
    #        dave2NO += tmpd2NO
    #        dave2IO += tmpd2IO

    #    dave2NO /= 31
    #    dave2IO /= 31

    #    dave_NO.append(dave1NO - dave2NO)
    #    dave_IO.append(dave1IO - dave2IO)

    #    print(mod[i], dave_NO[-1], dave_IO[-1])
    



    #### single mode distance :

    for i in range(1):
        dmin1NO = 1000000
        dmin1IO = 1000000
        for j in range(1):
            tmpd1NO = deltaCDF(cfd_arr_NO[i], cfd_arr_IO[j])
            tmpd1IO = deltaCDF(cfd_arr_IO[i], cfd_arr_NO[j])
            if tmpd1NO < dmin1NO :
                dmin1NO = tmpd1NO
            if tmpd1IO < dmin1IO:
                dmin1IO = tmpd1IO


    print(dmin1NO, dmin1IO)



    Nnu = 804
    count_arr = np.ones(Nnu)
    evtID, time = loadData("/junofs/users/miaoyu/supernova/production/Data/2kpc/data_mod82702_cha1_mh2_val.root")
    nn = 0
    dis_dist = []
    evt_time = []
    
    with open("mod82702IO.txt", "w") as ff:

        for i in time:
            if nn == Nnu:
                binCent0, binCont0 = [], []
                nn = 0
                N0 = cdf(Tmin, Tmax, evt_time, count_arr)
                for t in np.arange(Tmin, Tmax, bin_width):
                    N = cdf(Tmin, t, evt_time, count_arr)
                    binCent0.append((t-Tmin) / (Tmax - Tmin))
                    binCont0.append(N / N0)

                tmp = deltaCDF(cfd_arr_NO[0], binCont0)
                dis_dist.append(tmp)


                evt_time = []
                continue
            else:
                evt_time.append(i)
                nn += 1
    
    dis_dist = np.array(dis_dist)
    




    fig, ax = plt.subplots()
    ax.hist(dis_dist, bins=30, color="gray")
    ax.vlines(np.mean(dis_dist), 0, 100, color="red", lw=2.5)
    ax.vlines(np.mean(dis_dist) - np.std(dis_dist), 0, 100, color="red", lw=2.5, linestyle="--")
    ax.vlines(np.mean(dis_dist) + np.std(dis_dist), 0, 100, color="red", lw=2.5, linestyle="--")
    ax.vlines(0.087259, 0, 100, color="black", lw=2.5)
    ax.set_xlabel("distance", fontsize=14)
    ax.set_ylabel("count", fontsize=14)
    plt.show()




def calcCFD():
    path = "../production/PDFs/2kpc/evistSpec_mod"
    mod = [ 81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502,92503,92700,92701,92702,92703,94000,94001,94002,94003]
    ich = 1
    cfd_arr_NO, cfd_arr_IO = [], []

    fig, ax = plt.subplots()
    cid = 0
    for imod in mod:
        for imh in [1, 2]:

            filename = path + str(imod) + "_cha" + str(ich) + ".root"
            histname = "hET_mod"+str(imod)+"_cha"+str(ich)+"_mh"+str(imh)
            print(filename, histname)
            
            binCent, binCont = projection(filename, histname, ich)


            if imh == 1:   # NO
                binCent1, binCont1 = [], []
                N0 = cdf(Tmin, Tmax, binCent, binCont)
                for t in np.arange(Tmin, Tmax, bin_width):
                    N = cdf(Tmin, t, binCent, binCont)
                    binCent1.append((t-Tmin) / (Tmax - Tmin))
                    binCont1.append(N / N0)
                binCent1 = np.array(binCent1)
                binCont1 = np.array(binCont1)
                cfd_arr_NO.append(binCont1)
                #ax.plot(binCent1, binCont1/N0, lstyles[imh-1],  lw=2.0, color=colors[cid], label="mod%d_cha%d_mh%d"%(imod, ich, imh))
            if imh == 2:   # IO
                binCent2, binCont2 = [], []
                N0 = cdf(Tmin, Tmax, binCent, binCont)
                for t in np.arange(Tmin, Tmax, bin_width):
                    N = cdf(Tmin, t, binCent, binCont)
                    binCent2.append((t-Tmin) / (Tmax - Tmin))
                    binCont2.append(N / N0)
                binCent2 = np.array(binCent2)
                binCont2 = np.array(binCont2)
                cfd_arr_IO.append(binCont2)
                #ax.plot(binCent1, binCont1/N0, lstyles[imh-1],  lw=2.0, color=colors[cid], label="mod%d_cha%d_mh%d"%(imod, ich, imh))
        cid += 1
        print("****************************************"+str(imod)+"**********************************************")



    #with open("cfd_allMod.txt", "w") as f:
    #    for i, j in zip(mod, cfd_arr_NO):
    #        for k in range(len(j)):
    #            f.write("1"+" " + str(i) +" "+ str(binCent1[k]) + " " + str( j[k]) )
    #            f.write("\n")
    #    for i, j in zip(mod, cfd_arr_IO):
    #        for k in range(len(j)):
    #            f.write("2"+" " + str(i) +" "+ str(binCent1[k]) + " " + str( j[k]) )
    #            f.write("\n")

    Nnu = 804
    count_arr = np.ones(Nnu)
    evtID, time = loadData("/junofs/users/miaoyu/supernova/production/Data/2kpc/data_mod82702_cha1_mh2_val.root")
    nn = 0
    dis_dist = []
    evt_time = []
    for i in time:
        if nn == Nnu:
            binCent0, binCont0 = [], []
            nn = 0
            N0 = cdf(Tmin, Tmax, evt_time, count_arr)
            for t in np.arange(Tmin, Tmax, bin_width):
                N = cdf(Tmin, t, evt_time, count_arr)
                binCent0.append((t-Tmin) / (Tmax - Tmin))
                binCont0.append(N / N0)
            
            for mm in mod:
                tmp = deltaCDF(cfd_arr_NO[mm], binCont0)
                dis_dist.append(tmp)

            evt_time = []
            continue
        else:
            evt_time.append(i)
            nn += 1
    
    dis_dist = np.array(dis_dist)


if __name__ == "__main__" :
    calcCFD()
