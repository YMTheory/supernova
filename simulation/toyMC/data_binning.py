import numpy as np
import matplotlib.pyplot as plt
import uproot as up
import argparse


if __name__ == "__main__":
    
    MO      = "NO"
    model   = "Garching"
    modelNo = 82507
    Ethr    = 0.15
    fitTmin = -20
    fitTmax = 20
    fitTwidth = 0.01
    dist    = 10
    startevt = 0
    endevt  = 100000
    cha     = "IBD"
    
    parser = argparse.ArgumentParser(description='Arguments of SNNu analyser.')
    parser.add_argument('--model',      type=str,   default='Garching', help='Model name of SNNu.')
    parser.add_argument('--modelNo',    type=int,   default=82703,      help="modelNo")
    parser.add_argument('--MO',         type=str,   default="NO",       help="Mass ordering for the dataset.")
    parser.add_argument('--Ethr' ,      type=float, default=0.15,       help="Detection threshold for pES channel, unit MeV.")
    parser.add_argument('--fitTmin',    type=int,   default=-20,        help="Minimum fitting time with unit ms.")
    parser.add_argument('--fitTmax',    type=int,   default=20,         help="Maximum fitting time with unit ms.")
    parser.add_argument('--fitTwidth',  type=float, default=0.01,       help="binning width with unit ms.")
    parser.add_argument('--dist',       type=int,   default=10,         help="CCSNe distance")
    parser.add_argument('--channel',    type=str,   default="IBD",      help="Detection channel.")
    parser.add_argument("--start",      type=int,   default=0,          help="Start event for fit.")
    parser.add_argument("--end",        type=int,   default=0,          help="End event for fit.")
    args = parser.parse_args()

    model   = args.model
    modelNo = args.modelNo
    Ethr    = args.Ethr
    MO      = args.MO
    fitTmin = args.fitTmin
    fitTmax = args.fitTmax
    fitTwidth = 0.01
    dist    = args.dist
    startevt = args.start
    endevt  = args.end
    cha     = args.channel

    Nbins = int((fitTmax - fitTmin) / fitTwidth)

    datafile = f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha}_binneddata_{MO}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_merger.root"
    print(f"Loading datafile {datafile}...")
    with up.open(datafile) as f:
        nuTime = f["binned"]["TbinCont"].array(entry_start=startevt, entry_stop=endevt)
        nuNUM = []
        binConts = []
        for evt, nuTime_oneEvt in enumerate(nuTime):
            if evt % 100 == 0:
                print(f"**********Processing event {evt} ... **************")
            nuNUM.append(len(nuTime_oneEvt))
            
            conts, _ = np.histogram(nuTime_oneEvt, bins=Nbins, range=(fitTmin, fitTmax))
            binConts.append(conts)


    nuNUM = np.array(nuNUM)
    binConts = np.array(binConts)
    
    if 0:
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.hist(nuNUM, bins=100, color='blue', edgecolor="black")
        ax.set_xlabel(r"Number of $\nu$ in fitting time window", fontsize=14)
        ax.set_ylabel("A.U.", fontsize=14)
        plt.tight_layout()
        plt.show()
    
    outfile = f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha}_binneddata_{MO}_{dist}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning.root"
    with up.recreate(outfile) as f:
        f["binned"] = {"TbinConts" : binConts}
    
    
    



