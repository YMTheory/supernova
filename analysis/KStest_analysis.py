import uproot as up
import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__" :

    dist = int(sys.argv[1])

    testNO, testIO = [], []
    mod_arr = [81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502, 92503,92700,92701,92702,92703,94000,94001,94002,94003]
    for imod in mod_arr:
        if dist == 10000:
            filename = "/junofs/users/miaoyu/supernova/analysis/submit/KStest_mod%d_stat10000.root"%(imod)
        elif dist == 2:
            filename = "/junofs/users/miaoyu/supernova/analysis/submit/KStest_mod%d.root"%(imod)
        else:
            filename = "/junofs/users/miaoyu/supernova/analysis/submit/KStest_mod%d_%dkpc.root"%(imod, dist)
        ff = up.open(filename)
        tmpNO = ff["KS"]["NO"].array()
        tmpIO = ff["KS"]["IO"].array()

        testNO.extend(tmpNO)
        testIO.extend(tmpIO)

    testNO = np.array(testNO)
    testIO = np.array(testIO)

    print("Noraml Ordering Efficiency : ", np.sum(testNO>0) / float(len(testNO)))
    print("Inverted Ordering Efficiency: ", np.sum(testIO<0)/float(len(testIO)))

    fig, ax = plt.subplots()
    ax.hist(testNO, range=(-1, 1), bins=100, color="blue",    edgecolor="black", label="Noraml Ordering")
    ax.hist(testIO, range=(-1, 1), bins=100, color="crimson", edgecolor="black", label="Inverted Ordering")

    ax.vlines(0, 0, 800, color="black", linestyle="--")

    ax.set_xlabel("Test statistics", fontsize=13)
    ax.legend(prop={"size":13})


    plt.tight_layout()
    plt.savefig("KStest_%dkpc.pdf"%dist)
    plt.show()









