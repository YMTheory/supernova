#!/usr/bin/env python
# coding: utf-8


import fitter as fitter

from channel_analyser import channel

import numpy as np
import matplotlib.pyplot as plt
plt.style.use("science")

## Channel configuration:

model = "Garching"
modelNo = 82703
dist = 10

channels = {}
channels["eES"] = channel("eES", "IO", "Garching", modelNo, 0.20, fitTmin=-0.02, fitTmax=0.02, dist=dist, exp="JUNO")
channels["IBD"] = channel("IBD", "IO", "Garching", modelNo, 0.20, fitTmin=-0.02, fitTmax=0.02, dist=dist, exp="JUNO")
channels["pES"] = channel("pES", "IO", "Garching", modelNo, 0.10, fitTmin=-0.02, fitTmax=0.02, dist=dist, exp="JUNO")

for cha in channels.values():
    # set pdf file names:
    cha.setNOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_NO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")
    cha.setIOPdfFilePath(f"/junofs/users/miaoyu/supernova/simulation/C++/jobs/{model}{modelNo}_PDF_IO_10kpc_{cha.name}_{cha.Ethr:.2f}MeV_newshortPDF_v2.root")

    cha._load_pdf()


Garchings = [81123, 82503, 82703, 84003, 91123, 92503, 92703, 94003]
Garchings_model = [r"LS220, $11.2M_\odot$", r"LS220, $25M_\odot$", r"LS220, $27M_\odot$", r"LS220, $40M_\odot$", r"Shen, $11.2M_\odot$", r"Shen, $25M_\odot$", r"Shen, $27M_\odot$", r"Shen, $40M_\odot$",]
#dchi2_NO = np.zeros((8, 8))
#for i, dM in enumerate(Garchings):
#    for j, pM in enumerate(Garchings):
#        dTNO, chi2NO, dTIO, chi2IO, dchi2 = fitter.scanning_randomModels_chain(channels.values(), "NO", mode="fix", dataModel=dM, pdfModel=pM, group="Garching")
#        print(dM, pM, dTNO, chi2NO, dTIO, chi2IO, dchi2)
#        #dchi2_NO[i, j] = dchi2

#dchi2_IO = np.zeros((8, 8))
#for i, dM in enumerate(Garchings):
#    for j, pM in enumerate(Garchings):
#        dTNO, chi2NO, dTIO, chi2IO, dchi2 = fitter.scanning_randomModels_chain(channels.values(), "IO", mode="fix", dataModel=dM, pdfModel=pM, group="Garching")
#        print(dM, pM, dchi2)
#        print(dM, pM, dTNO, chi2NO, dTIO, chi2IO, dchi2)
#        dchi2_IO[i, j] = dchi2

if 1:
    print("Drawing...")
    dchi2_NO, dchi2_IO = np.zeros((8, 8)), np.zeros((8, 8))
    arrNO = np.loadtxt("models1_NO.txt")
    n = 0
    for i in range(8):
        for j in range(8):
            dchi2_NO[i, j] = arrNO[:, 2][n] - arrNO[:, 4][n]
            n+=1
    arrIO = np.loadtxt("models1_IO.txt")
    n = 0
    for i in range(8):
        for j in range(8):
            dchi2_IO[i, j] = arrIO[:, 2][n] - arrIO[:, 4][n]
            n+=1
    print(dchi2_NO)
    print(dchi2_IO)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    im0 = ax[0].imshow(dchi2_NO)
    cb0 = plt.colorbar(im0, ax=ax[0], shrink=0.6)
    ax[0].set_xticks(np.arange(0, 8, 1), Garchings_model, rotation=45)
    ax[0].set_yticks(np.arange(0, 8, 1), Garchings_model)
    ax[0].set_xlabel("PDF model", fontsize=18)
    ax[0].set_ylabel("Data model", fontsize=18)
    ax[0].text(8.1, 0, r"$\Delta t$", fontsize=18)
    #ax[0].text(8.1, 0, r"$\Delta\chi^2$", fontsize=18)
    ax[0].tick_params(axis="both", labelsize=14)
    ax[0].set_title("NO", fontsize=18)
    cb0.ax.tick_params(labelsize=14)
    
    im1 = ax[1].imshow(dchi2_IO)
    cb1 = plt.colorbar(im1, ax=ax[1], shrink=0.6)
    ax[1].set_xticks(np.arange(0, 8, 1), Garchings_model, rotation=45)
    ax[1].set_yticks(np.arange(0, 8, 1), Garchings_model)
    ax[1].set_xlabel("PDF model", fontsize=18)
    ax[1].set_ylabel("Data model", fontsize=18)
    ax[1].text(8.1, 0, r"$\Delta t$", fontsize=18)
    #ax[1].text(8.1, 0, r"$\Delta\chi^2$", fontsize=18)
    ax[1].tick_params(axis="both", labelsize=14)
    ax[1].set_title("IO", fontsize=18)
    cb1.ax.tick_params(labelsize=14)
    
    plt.tight_layout()
    #plt.savefig("test_Garchingmodels_uncertainty.pdf")
    plt.savefig("discrimination_time.pdf")
    plt.show()
