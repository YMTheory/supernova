from xsec.NuE_XS import NuE_XS

if __name__ == "__main__" :
    nue_xs = NuE_XS()
    print(nue_xs.diffXS(10, 4, 1))
