import numpy as np
import matplotlib.pyplot as plt

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

        
