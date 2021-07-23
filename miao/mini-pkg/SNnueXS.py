import numpy as np
#import matplotlib.pyplot as plt


def differentialXS(E, T, tp):
    me = 0.511
    if E<(T/2+np.sqrt(T*(T+me))/2.):
        return 0
    diffxs = 0
    sin2_thetaW = 0.23
    epsi_p, epsi_m = 0, 0
    indexG = 1.732e-44
    if tp == 2 or tp==4:
        epsi_p = -1*sin2_thetaW
        epsi_m = 0.5 - sin2_thetaW
        diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))
        return diffxs

    if tp==3 or tp==5:
        epsi_p = -1*sin2_thetaW
        epsi_m = 0.5-sin2_thetaW
        diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))
        return diffxs

    if tp == 0:
        epsi_p = -1*sin2_thetaW
        epsi_m = -0.5-sin2_thetaW
        diffxs = indexG*(epsi_m*epsi_m+epsi_p*epsi_p*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))
        return diffxs

    if tp ==1 :
        epsi_p = -1*sin2_thetaW
        epsi_m = -0.5-sin2_thetaW
        diffxs = indexG*(epsi_p*epsi_p+epsi_m*epsi_m*np.power((1-T/E),2)-epsi_p*epsi_m*me*T/(E*E))

        return diffxs

    return 0 



def totalXS(E, tp):
    totalxs = 0
    sin2_thetaW = 0.23
    indexG = 1.732e-44
    if tp == 2 or tp == 4:
        totalxs = (indexG/4.*E) * (1-4*sin2_thetaW+16/3*np.power(sin2_thetaW, 2))
        return totalxs
    if tp ==3 or tp == 5:
        totalxs = (indexG/4.*E) * (1/3. - 4/3.*sin2_thetaW+16./3.*np.power(sin2_thetaW, 2))
        return totalxs
    if tp == 0:
        totalxs = (indexG/4.*E) * (1+4*sin2_thetaW+16./3.*np.power(sin2_thetaW, 2))
        return totalxs
    if tp == 1:
        totalxs = (indexG/4.*E) * (1/3. + 4/3.*sin2_thetaW+16./3.*np.power(sin2_thetaW, 2))
        return totalxs
    return 0.



"""
def xsPlot():

    E = np.arange(0.1, 20, 0.1)
    ty = [0, 1, 2, 3]
    tyname = [r"$\nu_e$", r"$\nu_{\bar{e}}$", r"$\nu_x$", r"$\nu_{\bar{x}}$"]
    xsec = [[]  for i in range(4)]
    cn = 0
    for t in ty:
        for i in E:
            xsec[cn].append(totalXS(i, t))
        plt.plot(E,  xsec[cn], label=tyname[cn])
        cn += 1

    plt.xlabel(r"$E_\nu$/MeV")
    plt.ylabel("totalXS")
    plt.grid(True)
    plt.legend()
    plt.savefig("nue_totalxs.pdf")


if __name__ == "__main__":
    xsPlot()

"""
        
