import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
sys.path.append("/junofs/users/miaoyu/supernova/simulation/toyMC/script/")
import plotConfig

def readCSV(filename, start, end):
    df1 = pd.read_csv(filename)
    Tbest = df1["Tbest"][start:end]
    locMin = df1["locMin"][start:end]
    return Tbest, locMin



def load_profiles(MO, start=0, end=100, param="chi2"):

    Nevt = end - start
    masses = np.arange(0.0, 2.0, 0.1)
    if param == "chi2":
        chi2_array = np.zeros((Nevt, len(masses)))
    elif param == "deltat":
        deltat_array = np.zeros((Nevt, len(masses)))


    csvfilepath = os.getenv("CSVFILEPATH")
    for i, mass in enumerate(masses):
        filename = f"{csvfilepath}Garching82703_10.0kpc_20.0kton_{MO}_pESFalseeESTrueIBDTrue_nuMass{mass:.1f}eV_0.10MeV_fitTmin-0.020sfitTmax0.020s_noC14_tot_unbinnedData_unbinnedNLL_PosiToyDataTobs2D_JUNO_limitdTtoAvoidZeroProb.csv"
        tmp_Tbest, tmp_locMin = readCSV(filename, start, end)
        for j in range(Nevt):
            if param == "chi2":
                chi2_array[j, i] = 2*tmp_locMin[j]
            elif param == "deltat":
                deltat_array[j, i] = tmp_Tbest[j]

    if param == "chi2":
        return chi2_array
    elif param == "deltat":
        return deltat_array


def linear_interpolate_crossing(dm, chi2, critical_value=2.706):
    for i in range(len(dm)-1):
        if chi2[i] <= critical_value and chi2[i+1] >= critical_value:
            ## suppose this is a rising crossing point:
            dy0 = chi2[i+1] - chi2[i]
            dx0 = dm[i+1] - dm[i]
            dy = critical_value - chi2[i]
            dx = dx0 * dy / dy0
            cross = dm[i] + dx
            return cross
    else:
        return -100


def discard_abnormalFit(arr):
    if np.min(arr) < -9000:
        return True
    else:
        return False


def median_sensitivity(start=0, end=100000):
     
    nuMass = np.arange(0, 2.0, 0.1)
    chi2_profilesNO = load_profiles("NO", start=start, end=end, param="chi2")
    chi2_profilesIO = load_profiles("IO", start=start, end=end, param="chi2")

    ievt, Ndiscaded_NO, Ndiscaded_IO = 0, 0, 0
    crossing_points_NO, crossing_points_IO = [], []
    for chi2NO, chi2IO in zip(chi2_profilesNO, chi2_profilesIO):
        if discard_abnormalFit(chi2NO):
            Ndiscaded_NO += 1
        else :
            chi2NO = chi2NO - np.min(chi2NO)
            crossNO = linear_interpolate_crossing(nuMass, chi2NO)
            crossing_points_NO.append(crossNO)

        if discard_abnormalFit(chi2IO):
            Ndiscaded_IO += 1
        else:
            chi2IO = chi2IO - np.min(chi2IO)
            crossIO = linear_interpolate_crossing(nuMass, chi2IO)
            crossing_points_IO.append(crossIO)
        
        ievt += 1
        if ievt % 100 == 0:
            print(f"Processing event ID = {ievt} => Discarding {Ndiscaded_NO} NO events and {Ndiscaded_IO} IO events.")

    # mask abnoraml events:
    crossing_points_NO = np.array(crossing_points_NO)
    crossing_points_IO = np.array(crossing_points_IO)
    masked_crossing_points_NO = ma.masked_where(crossing_points_NO == -100, crossing_points_NO)
    masked_crossing_points_IO = ma.masked_where(crossing_points_IO == -100, crossing_points_IO)

    plot = True
    if plot:
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))

        ax[0].hist(masked_crossing_points_NO, bins=50, range=(0, 3), color="blue", edgecolor="black")
        ax[0].vlines(ma.median(masked_crossing_points_NO), 0, 4000, linestyle="--", color="white", linewidth=2)
        ax[0].text(1.9, 5000, f"{ma.median(masked_crossing_points_NO):.2f}"+r"$\pm$"+f"{ma.std(masked_crossing_points_NO):.2f}eV", fontsize=16) 
        ax[1].hist(masked_crossing_points_IO, bins=50, range=(0, 3), color="blue", edgecolor="black")
        ax[1].vlines(ma.median(masked_crossing_points_IO), 0, 4000, linestyle="--", color="white", linewidth=2)
        ax[1].text(1.9, 5000, f"{ma.median(masked_crossing_points_IO):.2f}"+r"$\pm$"+f"{ma.std(masked_crossing_points_IO):.2f}eV", fontsize=16) 

        plotConfig.setaxis(ax[0], xlabel=r"$m_\nu$ [eV]", ylabel="", yticksize=0, title="NO")
        plotConfig.setaxis(ax[1], xlabel=r"$m_\nu$ [eV]", ylabel="", yticksize=0, title="IO")

        plt.tight_layout()
        plt.savefig("./plots/medianSens_Garching82703_10.0kpc_nopES_limitdTtoAvoidZeroProb.pdf")
        plt.show()



if __name__ == "__main__":
    median_sensitivity()












