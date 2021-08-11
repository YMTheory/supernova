pdfPathPrefix = "./etSpec/fineSpec/"
toyDataPathPrefix = "./dataset/fineSpec/"
fitResPathPrefix = "./dataset/"
fitPlotPathPrefix = "./spectra/fineSpec/"
LogPathPrefix = "./log/"
ShPathPrefix = "./job/"

chaName = ["nup", "nue", "IBD"]

def pdfFileName(modelNum, dist, chaname, fitNuMass, fitMH, nuType, Evismax, no):
    ## fitting pdf filename
    MO = ["noOsc", "NO", "IO"]
    prefix = pdfPathPrefix + "%dkpc/" %(dist)
    prefix += MO[int(fitMH)] + "/"
    filename = prefix + 'TEvisPDF_mod%d_cha%d%s_nuType%d_mh%d'%(modelNum, chaname, chaName[chaname], nuType, fitMH)
    filename = filename + '_mNu%2.1feV_%2.1fkpc_0.1s_Evismax%.2f_No%d.root'%(fitNuMass, dist, Evismax, no)
    print("-----> PDF filename: %s" %filename)

    return filename


def pdfSumFileName(modelNum, dist, chaname, fitNuMass, fitMH, nuType, Evismax):
    ## fitting pdf filename
    MO = ["noOsc", "NO", "IO"]
    prefix = pdfPathPrefix + "%dkpc/" %(dist)
    prefix += MO[int(fitMH)] + "/"
    filename = prefix + 'TEvisPDF_mod%d_cha%d%s_nuType%d_mh%d'%(modelNum, chaname, chaName[chaname], nuType, fitMH)
    filename = filename + '_mNu%2.1feV_%2.1fkpc_0.1s_Evismax%.2f_Sum.root'%(fitNuMass, dist, Evismax)
    print("-----> PDF filename: %s" %filename)

    return filename





def pdfEvisHist2DName(modelNum, chaname, fitMH):
    histname = "hET_mod%d_cha%d_mh%d" %(modelNum, chaname, fitMH)
    return histname


def pdfEvHist2DName(modelNum, chaname, fitMH):
    histname = "hENuT_mod%d_cha%d_mh%d" %(modelNum, chaname, fitMH)
    return histname


def pdfBkgHist2DName(modelNum, chaname, fitMH):
    histname = "hETBKG_mod%d_cha%d_mh%d" %(modelNum, chaname, fitMH)
    return histname

def pdfGenLogName(modelNum, chaname, nuType, nuMass, MH, dist, no):
    logname = LogPathPrefix + "log-genPDFmiao-mod%d-cha%d-nuType%d-nuMass%2.1feV-MH%d-%2.1fdist-No%d.txt"%(modelNum, chaname, nuType, nuMass, MH, dist, no)
    return logname

def pdfGenShName(modelNum, chaname, nuType, nuMass, MH, dist, no):
    logname = ShPathPrefix + "run-genPDFmiao-mod%d-cha%d-nuType%d-nuMass%2.1feV-MH%d-%2.1fdist-No%d.sh"%(modelNum, chaname, nuType, nuMass, MH, dist, no)
    return logname


def dataFileName(modelNum, dist, chaname, dataNuMass, dataMH, Ethr, group):
    filename = toyDataPathPrefix + "%2.1fkpc/TEvisDATA_" %(dist)
    filename = filename + 'mod%d_cha%d%s_mh%d' %(modelNum, chaname, chaName[chaname], dataMH)
    filename = filename + '_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25_Ethr%2.1fMeV_group%d.root' %(dataNuMass, dist, Ethr, group)

    print("-----> toyData filename: %s" %filename)
    return filename
    
    

def fitResFileName(modelNum, chaname, dataMH, dataNuMass, fitMH, fitNuMass, dist, Ethr, group):
    filename = fitResPathPrefix + "%2.1fkpc/TEvisDATA_" %(dist)
    filename = filename + "mod%d_cha%d%s_dataMH%d_datamNu%2.1feV_fitMH%d_fitmNu%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_1DFit" %(modelNum, chaname, chaName[chaname], dataMH, dataNuMass, fitMH, fitNuMass, dist, Ethr, group)

    rawfile = filename + "Raw.txt"
    sumfile = filename + "Summary.txt"

    return rawfile, sumfile



def fitResPdfName(modelNum, chaname, dataMH, dataNuMass, fitMH, fitNuMass, dist, Ethr, group):
    filename = fitResPathPrefix + "%2.1fkpc/TEvisDATA_" %(dist)
    filename = filename + "mod%d_cha%d%s_dataMH%d_datamNu%2.1feV_fitMH%d_fitmNu%2.1feV_%2.1fkpc_Ethr%2.1fMeV_group%d_1DFit.pdf" %(modelNum, chaname, chaName[chaname], dataMH, dataNuMass, fitMH, fitNuMass, dist, Ethr, group)
    return filename

