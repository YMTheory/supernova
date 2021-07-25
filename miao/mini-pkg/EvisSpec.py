import numpy as np
import matplotlib.pyplot as plt
import SNGarchingIntegFcn as gar
from SNnumGarchingSrc import SNnumGarchingSrc
import SNnueXS as nuexs
from detector import detector
import boost_histogram as bh
from EvEneSpec import EvEneSpec
import ROOT

libSNsimDir = '/junofs/users/miaoyu/supernova/wenlj/simulation/lib/libSNsim.so'
#libSNsimDir = os.getenv('SNCODEDIR') + '/simulation/lib/libSNsim.so'
ROOT.gSystem.Load(libSNsimDir)
print('Load', libSNsimDir)


snDet = ROOT.SNdetect.instance()
modelNum = 82503
snDet.setSrcModel(modelNum)
modelSrc = snDet.getPointerSrc()

# -- configurations -- #
snDet.initChannel(1)
modelSrc.setSNDistance(10)
Ethr = 0.2
snDet.getPointerEffectLS().setThresholdE(Ethr)
snDet.initFCN()


print("------> eES channel w/o osc total event number = %d" %snDet.getEventAboveEthrVis(0.2, -1, 0))
print("------> eES channel w/  NO  total event number = %d" %snDet.getEventAboveEthrVis(0.2, -1, 1))
print("------> eES channel w/  IO  total event number = %d" %snDet.getEventAboveEthrVis(0.2, -1, 2))






"""
# for eES channel only ...
Evismin, Evismax, step_Evis = 0.2, 50, 0.05
Evmin,   Evmax,   step_Ev   = 0.0, 60, 0.2
nbins_Ev   = round( (Evmax-Evmin)/step_Ev )
nbins_Evis = round( (Evismax-Evismin)/step_Evis )

hist_Eve = bh.Histogram(bh.axis.Regular(bins=nbins_Ev, start=Evmin, stop=Evmax))
hist_Evebar = bh.Histogram(bh.axis.Regular(bins=nbins_Ev, start=Evmin, stop=Evmax))
hist_Evx = bh.Histogram(bh.axis.Regular(bins=nbins_Ev, start=Evmin, stop=Evmax))
hist_Evis = bh.Histogram(bh.axis.Regular(bins=nbins_Evis, start=Evismin, stop=Evismax))
hist_Evis_NO = bh.Histogram(bh.axis.Regular(bins=nbins_Evis, start=Evismin, stop=Evismax))
hist_Evis_IO = bh.Histogram(bh.axis.Regular(bins=nbins_Evis, start=Evismin, stop=Evismax))
x_Ev = hist_Eve.axes.centers[0]
x_Evis = hist_Evis.axes.centers[0]

det = detector()
Npar = det.getNumberOfProton() + 6*det.getNumberOfCarbon()
spec = EvEneSpec(1, 12, 1)

#nuRate = 0
## neutrino rate :
#for i in range(len(x_Ev)):
#    EvTmp = x_Ev[i]
#    fluence = spec.KRJ_dFdE(EvTmp)
#    txs = nuexs.totalXS(EvTmp, 0)
#    nevt = fluence * txs * Npar
#    print(EvTmp, fluence, txs)
#    nuRate += nevt * step_Ev
#
#print("Total Ev rate = %d"%nuRate)




# neutrino energy spectra :
# nu_e
cont_Evis = [ 0 for i in range(nbins_Evis)]
cont_Evis_miao = [ 0 for i in range(nbins_Evis)]
cont_Evis_NO = [ 0 for i in range(nbins_Evis)]
cont_Evis_IO = [ 0 for i in range(nbins_Evis)]
SNGar = SNnumGarchingSrc(82503, 10)
for i in range(len(x_Ev)):
    EvTmp = x_Ev[i]
    for j in range(len(x_Evis)):
        EvisTmp = x_Evis[j]
        fluence = SNGar.oneSNFluenceDetMH(EvTmp, 0, 0)
        fluence0 = spec.KRJ_dFdE(EvTmp)
        fluence1 = SNGar.oneSNFluenceDetMH(EvTmp, 0, 1)
        fluence2 = SNGar.oneSNFluenceDetMH(EvTmp, 0, 2)
        hist_Eve[i] = fluence
        dxs = nuexs.differentialXS(EvTmp, EvisTmp, 0)
        if fluence > 0 and dxs > 0:
            cont_Evis[j] += fluence * dxs * Npar
            cont_Evis_NO[j] += fluence1 * dxs * Npar
            cont_Evis_IO[j] += fluence2 * dxs * Npar
            cont_Evis_miao[j] += fluence0 * dxs * Npar

# nu_e
for i in range(len(x_Ev)):
    EvTmp = x_Ev[i]
    for j in range(len(x_Evis)):
        EvisTmp = x_Evis[j]
        fluence = SNGar.oneSNFluenceDetMH(EvTmp, 1, 0)
        fluence1 = SNGar.oneSNFluenceDetMH(EvTmp, 1, 1)
        fluence2 = SNGar.oneSNFluenceDetMH(EvTmp, 1, 2)
        hist_Evebar[i] = fluence
        dxs = nuexs.differentialXS(EvTmp, EvisTmp, 1)
        if fluence > 0 and dxs > 0:
            #cont_Evis[j] += fluence * dxs * Npar
            cont_Evis_NO[j] += fluence1 * dxs * Npar
            cont_Evis_IO[j] += fluence2 * dxs * Npar

## nu_e
for i in range(len(x_Ev)):
    EvTmp = x_Ev[i]
    for j in range(len(x_Evis)):
        EvisTmp = x_Evis[j]
        fluence = SNGar.oneSNFluenceDetMH(EvTmp, 2, 0)
        fluence1 = SNGar.oneSNFluenceDetMH(EvTmp, 2, 1)
        fluence2 = SNGar.oneSNFluenceDetMH(EvTmp, 2, 2)
        hist_Evx[i] = fluence
        dxs = nuexs.differentialXS(EvTmp, EvisTmp, 2)
        if fluence > 0 and dxs > 0:
            #cont_Evis[j] += fluence * dxs * Npar
            cont_Evis_NO[j] += 4 * fluence1 * dxs * Npar
            cont_Evis_IO[j] += 4 * fluence2 * dxs * Npar

for i in range(nbins_Evis):
    hist_Evis[i] = cont_Evis[i]
    hist_Evis_NO[i] = cont_Evis_NO[i]
    hist_Evis_IO[i] = cont_Evis_IO[i]

cont_Evis = np.array(cont_Evis)
cont_Evis_NO = np.array(cont_Evis_NO)
cont_Evis_IO = np.array(cont_Evis_IO)
cont_Evis_miao = np.array(cont_Evis_miao)



print("Total Event Entries w/o oscillation mine = %d" %(np.sum(cont_Evis_miao)*step_Evis))
print("Total Event Entries w/o oscillation = %d" %(np.sum(cont_Evis)*step_Evis))
print("Total Event Entries w/ oscillation NO = %d" %(np.sum(cont_Evis_NO)*step_Evis))
print("Total Event Entries w/ oscillation IO = %d" %(np.sum(cont_Evis_IO)*step_Evis))




# visualisation : 
plt.figure(0)
plt.plot(x_Ev, hist_Eve.values(), "-", label=r"$\nu_e$")
plt.plot(x_Ev, hist_Evebar.values(), "-", label=r"$\bar\nu_{e}$")
plt.plot(x_Ev, hist_Evx.values(), "-", label=r"$\nu_x$")
plt.legend()
plt.xlabel(r"$E_\nu$/MeV")
plt.ylabel(r"dN/dE [MeV$^{-1}$]")
plt.savefig("./outputs/Energy_spectra_allType_noOsc.pdf")

#plt.figure(1)
#plt.plot(x_Evis, x_Evis*hist_Evis.values(), "-")
#plt.xlabel(r"$E_d$/MeV")
#plt.ylabel(r"EdN/dE [MeV$^{-1}$]")
#plt.ylim(0.1, 100000)
#plt.semilogx()
#plt.semilogy()
#plt.savefig("eES_EddNoverdE.pdf")



plt.figure(1)
plt.plot(x_Evis, x_Evis*hist_Evis.values(),    "-",  color="blue", label="w/o osc")
plt.plot(x_Evis, x_Evis*cont_Evis_miao,    "-",  color="red", label="w/o osc")
plt.plot(x_Evis, x_Evis*hist_Evis_NO.values(), "--", color="blue", label="w/ osc, NO")
plt.plot(x_Evis, x_Evis*hist_Evis_IO.values(), "-.", color="blue", label="w/ osc, IO")
plt.xlabel(r"$E_d$/MeV")
plt.ylabel(r"EdN/dE [MeV$^{-1}$]")
#plt.ylim(0.1, 100000)
plt.semilogx()
plt.semilogy()
plt.legend()
#plt.savefig("./outputs/eES_EddNoverdE.pdf")
plt.show()
"""
