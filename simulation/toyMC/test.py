from channel_analyser import channel

pES = channel("pES", "NO", "Garching", 82703, 0.15)
pES.setNevtPerFile(100000)
pES.setStartEvtId(2000)
pES.setEndEvtId(3000)
pES.setDataFilePath(f"/junofs/users/miaoyu/supernova/simulation/toyMC/scale10/Garching82703_pES_data_NO_10kpc_thr0.15MeV_Tmin10msTmax50ms_merger.root")
pES._load_data()

