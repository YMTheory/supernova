from channel_analyser import channel

if __name__ == "__main__":

    MO          = "NO"
    model       = "Burrows"
    modelNo     = 12
    Ethr        = 0.10
    fitTmin     = 10
    fitTmax     = 45
    fileNo      = 0
    dist        = 10

    channels    = {}
    channels["pES"] = channel("pES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax = fitTmax, fileNo=fileNo, dist=dist)
    channels["IBD"] = channel("IBD", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax = fitTmax, fileNo=fileNo, dist=dist)
    channels["eES"] = channel("eES", MO, model, modelNo, Ethr, fitTmin=fitTmin, fitTmax = fitTmax, fileNo=fileNo, dist=dist)

    for cha in channels.values():
        cha.setDataFilePath(f"/afs/ihep.ac.cn/users/m/miaoyu/junofs/supernova/simulation/toyMC/scale1_poisson/{model}{modelNo}_{cha.name}_binneddata_{MO}_{dist:.1f}kpc_thr{Ethr:.2f}MeV_Tmin{fitTmin}msTmax{fitTmax}ms_binning.root")
        cha._load_data_ak()
        print(cha.data_array)
