import numpy as np

def read_data_fromHL(filename):
    cc = 0
    T, N = [], []
    with open(filename) as f:
        nn = 0
        for lines in f.readlines():
            if cc == 6:
                T.append(t)
                N.append(nn)
                cc = 0
                nn = 0
            line = lines.strip("\n")
            data = line.split(" ")
            nn += float(data[2])
            t = float(data[0])
            cc += 1

    T = np.array(T) * 1000
    N = np.array(N) / 1000.
    return T, N


import ROOT

if __name__ == "__main__" :
    #for imod in [81120,81121,81122,81123,82500,82501,82502,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502,92503,92700,92701,92702,92703,94000,94001,94002,94003]:
    for imod in [82503]:
        for imh in [1, 2]:
            for cha in [ 0, 1, 2]:
            #for cha in [0]:
                if cha == 1 or cha == 2:
                    filename = "/junofs/users/miaoyu/supernova/wenlj/simulation/examples/submit/Gar82503_cha%d_MH%d_Ethr%.2fMeV.txt"%(cha, imh, 0.20)
                    print(filename)
        
                    T, N = read_data_fromHL(filename)
                    #h1 = ROOT.TH1D("h1", "visible energy spectrum rate", 60, -0.5, 60.5)
                    h1 = ROOT.TH1D("h1", "visible energy spectrum rate", 60, -30.5, 30.5)

                    for i in range(len(T)):
                        h1.SetBinContent(i+1, N[i])

                    if imh == 1:
                        MO = "NO"
                    elif imh == 2:
                        MO = "IO"
                    chaname = ["pES", "eES", "IBD"]
                    f = ROOT.TFile("Garching82503_PDF_%s_10kpc_%s_%.2fMeV.root"%( MO, chaname[cha], 0.20), "recreate")
                    h1.Write()
                    f.Close()

                    del h1
                    del f


                if cha == 0:
                    for Ethr in [0.10, 0.15, 0.20]:
                        filename = "/junofs/users/miaoyu/supernova/wenlj/simulation/examples/submit/Gar82503_cha%d_MH%d_Ethr%.2fMeV.txt"%( cha, imh, Ethr)
                        print(filename)
        
                        T, N = read_data_fromHL(filename)
                        #h1 = ROOT.TH1D("h1", "visible energy spectrum rate", 60, -0.5, 60.5)
                        h1 = ROOT.TH1D("h1", "visible energy spectrum rate", 60, -30.5, 30.5)

                        for i in range(len(T)):
                            h1.SetBinContent(i+1, N[i])

                        if imh == 1:
                            MO = "NO"
                        elif imh == 2:
                            MO = "IO"
                        chaname = ["pES", "eES", "IBD"]
                        f = ROOT.TFile("Garching82503_PDF_%s_10kpc_%s_%.2fMeV.root"%(MO, chaname[cha], Ethr), "recreate")
                        h1.Write()
                        f.Close()

                        del h1
                        del f
    













