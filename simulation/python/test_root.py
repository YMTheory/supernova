from integration import integration

from tqdm import tqdm
import argparse
import numpy as np

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Arguments of SNNu generator")
    parser.add_argument('--channel', type=str,   default="NuP", help="detection channel name.")
    parser.add_argument('--tmin',    type=float, default=-0.03, help="start time in sec.")
    parser.add_argument('--tmax',    type=float, default=0.03,  help="end time in sec.")
    parser.add_argument('--Emin',    type=float, default=0.1,   help="minimum visible energy in MeV.")
    parser.add_argument('--Emax',    type=float, default=10.,   help="maximum visible energy in MeV.")
    args = parser.parse_args()

    cha = args.channel
    tmin = args.tmin
    tmax = args.tmax
    Emin = args.Emin
    Emax = args.Emax

    spec = integration(cha, SNmass=25, EoS="LS220")
    spec.build_fCore()
    print(f"=====> Integration time interval [{tmin}s, {tmax}s] and energy interval [{Emin}MeV, {Emax}MeV] *********")
    

    tarr = np.arange(tmin, tmax, 0.001)
    Evismin, Evismax = Emin, Emax

    arr = np.zeros((len(tarr), 2))

    for i in tqdm(range(len(tarr))):
        t = tarr[i]
        arr[i, 0] = t
        arr[i, 1] = spec.integrateEvisEvAtT(t, Evismin, Evismax, 0)


    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(tarr*1000, arr[:, 1]/1000, lw=2, label="NuE")
    ax.set_xlabel("post-bounce time [ms]", fontsize=14)
    ax.set_ylabel("Counts per ms", fontsize=14)
    #ax.set_ylabel(r"E$_\mathrm{vis}$ [MeV]", fontsize=14)
    ax.tick_params(axis="both", labelsize=13)
    ax.legend(prop={"size":12})
    #cb.set_label(r"s$^{-1}\cdot$MeV$^{-1}$", fontsize=14)
    plt.tight_layout()
    plt.show()
    





