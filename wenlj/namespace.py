pdfPathPrefix = "./etSpec/fineSpec/"
toyDataPathPrefix = "./dataset/fineSpec/"
fitResPathPrefix = "./dataset/"
fitPlotPathPrefix = "./spectra/fineSpec/"

chaName = ["nup", "nue", "IBD"]

def pdfFileName(modelNum, dist, chaname, fitNuMass, fitMH):
    ## fitting pdf filename
    filename = pdfPathPrefix + 'TEvisPDF_mod%d_cha%d%s_mh%d'%(modelNum, chaname, chaName[chaname], fitMH)
    filename = filename + '_mNu%2.1feV_%2.1fkpc_0.1s_Evmax25.root'%(fitNuMass, dist)
    print("-----> PDF filename: %s" %filename)

    return filename


def pdfSigHist2DName(modelNum, chaname, fitMH):
    histname = "hET_mod%d_cha%d_mh%d" %(modelNum, chaname, fitMH)
    return histname


def pdfBkgHist2DName(modelNum, chaname, fitMH):
    histname = "hETBKG_mod%d_cha%d_mh%d" %(modelNum, chaname, fitMH)
    return histname



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

