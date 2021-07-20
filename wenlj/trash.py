#     # ################################################
#     # ###    Draw the PDF
#     can = ROOT.TCanvas("2d fit", "2d fit", 1200, 400)
#     can.Divide(2)
#     can.cd(1)
#     can.Print(prefix+'_data%2.1feV_localMinNLL.pdf['%(nuMass1))

#     xframe = nuTime.frame(ROOT.RooFit.Title("test1"), ROOT.RooFit.Range('fullRange'),
#                         ROOT.RooFit.Bins(100))
#     xframe.GetYaxis().SetTitleOffset(1.4)
#     # dataPdf.plotOn(xframe, ROOT.RooFit.Range('fitRange'),
#     #                     ROOT.RooFit.DrawOption("F"),
#     #                     ROOT.RooFit.FillColor(ROOT.kOrange),
#     #                     ROOT.RooFit.MoveToBack() )
#     dataPdf.plotOn(xframe, ROOT.RooFit.Range('fullRange'),
#             ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#             ROOT.RooFit.LineColor(2))
#     xframe.Draw()

#     can.cd(2)
#     yframe = nuEnergy.frame(ROOT.RooFit.Title("test1"), ROOT.RooFit.Range('fullRange'),
#                         ROOT.RooFit.Bins(100))
#     yframe.GetYaxis().SetTitleOffset(1.4)
#     # fitdata.plotOn(yframe)
#     # fitPdf3.plotOn(xframe, ROOT.RooFit.ProjectionRange('fitRange'), 
#     #         ROOT.RooFit.Normalization(fitNevts, ROOT.RooAbsReal.NumEvent),
#     #         ROOT.RooFit.LineColor(2))
#     # fitPdf.plotOn(xframe, ROOT.RooFit.ProjectionRange('fitRange'), 
#     #         ROOT.RooFit.Normalization(fitNevts, ROOT.RooAbsReal.NumEvent),
#     #         ROOT.RooFit.LineColor(1))
#     # dataPdf.plotOn(yframe2, ROOT.RooFit.Range('fitRange'),
#     #                     ROOT.RooFit.DrawOption("F"),
#     #                     ROOT.RooFit.FillColor(ROOT.kOrange),
#     #                     ROOT.RooFit.MoveToBack() )
#     dataPdf.plotOn(yframe, ROOT.RooFit.Range('fullRange'),
#             ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#             ROOT.RooFit.LineColor(2))
#     yframe.Draw()

#     can.Print(prefix+'_data%2.1feV_localMinNLL.pdf'%(nuMass1))
#     ############################################
#     can.cd(1)
#     xframe = nuTime.frame(ROOT.RooFit.Title("test2"), ROOT.RooFit.Range('fullRange'),
#                     ROOT.RooFit.Bins(100))
#     dataPdf.plotOn(xframe, ROOT.RooFit.Range('fitRange'),
#         ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#         ROOT.RooFit.LineColor(2))
#     xframe.Draw()

#     can.cd(2)
#     yframe = nuEnergy.frame(ROOT.RooFit.Title("test2"), ROOT.RooFit.Range('fullRange'),
#                 ROOT.RooFit.Bins(100))
#     dataPdf.plotOn(yframe, ROOT.RooFit.Range('fitRange'),
#         ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#         ROOT.RooFit.LineColor(2))
#     yframe.Draw()

#     can.Print(prefix+'_data%2.1feV_localMinNLL.pdf'%(nuMass1))
#     ############################################
#     can.cd(1)
#     xframe = nuTime.frame(ROOT.RooFit.Title("test3"), ROOT.RooFit.Range('fitRange'),
#                     ROOT.RooFit.Bins(100))
#     dataPdf.plotOn(xframe, ROOT.RooFit.Range('fitRange'),
#         #ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#         ROOT.RooFit.LineColor(2))
#     xframe.Draw()

#     can.cd(2)
#     yframe = nuEnergy.frame(ROOT.RooFit.Title("test3"), ROOT.RooFit.Range('fitRange'),
#                     ROOT.RooFit.Bins(100))
#     dataPdf.plotOn(yframe, ROOT.RooFit.Range('fitRange'),
#         #ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#         ROOT.RooFit.LineColor(2))
#     yframe.Draw()

#     can.Print(prefix+'_data%2.1feV_localMinNLL.pdf'%(nuMass1))
#     ############################################
#     newArgSet = ROOT.RooArgSet(nuTime, nuEnergy)
#     integralValue = dataPdf.createIntegral(newArgSet,
#                             ROOT.RooFit.NormSet(newArgSet),
#                             ROOT.RooFit.Range('fitRange'))
#     print('integralValue', integralValue.getVal())

#     histTest = dataPdf.createHistogram("histTest", nuTime, 
#                 ROOT.RooFit.Binning(50), ROOT.RooFit.Range('fitRange'))
#     print('histTest.Integral()', histTest.Integral())

#     histTest2 = dataPdf.createHistogram("histTest2", nuEnergy, 
#                 ROOT.RooFit.Binning(50), ROOT.RooFit.Range('fitRange'))
#     print('histTest2.Integral()', histTest2.Integral())

#     can.cd(1)
#     xframe = nuTime.frame(ROOT.RooFit.Title("test4"), 
#                     ROOT.RooFit.Range('fitRange'),
#                     ROOT.RooFit.Bins(50))
#     xframe.GetYaxis().SetTitleOffset(1.4)
#     fitdata.plotOn(xframe)
#     # dataPdf.plotOn(xframe, ROOT.RooFit.Range('fitRange'),
#     #        #ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#     #        ROOT.RooFit.LineColor(2))
#     xframe.Draw()
#     histTest.Scale(nFitEvts)
#     histTest.Draw("same")

#     can.cd(2)
#     yframe = nuEnergy.frame(ROOT.RooFit.Title("test4"), 
#                         #ROOT.RooFit.Range('fitRange'),
#                         ROOT.RooFit.Bins(50))
#     yframe.GetYaxis().SetTitleOffset(1.4)
#     fitdata.plotOn(yframe)
#     dataPdf.plotOn(yframe, #ROOT.RooFit.Range('fitRange'),
#             ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#             ROOT.RooFit.LineColor(2))
#     yframe.Draw()
#     histTest2.Scale(nFitEvts)
#     histTest2.Draw("same")

#     can.Print(prefix+'_data%2.1feV_localMinNLL.pdf'%(nuMass1))

#     can.cd(1)
#     xframe = nuTime.frame(ROOT.RooFit.Title("Low statistics histogram pdf"), 
#                             ROOT.RooFit.Range('fitRange'),
#                             ROOT.RooFit.Bins(100) )
#     xframe.GetYaxis().SetTitleOffset(1.4)
#     fitdata.plotOn(xframe)
#     dataPdf.plotOn(xframe, ROOT.RooFit.NormRange('fitRange'),
#             #ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#             ROOT.RooFit.LineColor(2))
# #    fitPdf.plotOn(xframe, ROOT.RooFit.LineColor(2))
#     xframe.Draw()

#     can.cd(2)
#     yframe = nuEnergy.frame(ROOT.RooFit.Title("Low statistics histogram pdf"), 
#                             ROOT.RooFit.Range('fitRange'),
#                             ROOT.RooFit.Bins(100) )
#     yframe.GetYaxis().SetTitleOffset(1.4)
#     fitdata.plotOn(yframe)
#     dataPdf.plotOn(yframe, ROOT.RooFit.NormRange('fitRange'),
#             #ROOT.RooFit.Normalization(nFitEvts, ROOT.RooAbsReal.NumEvent),
#             ROOT.RooFit.LineColor(2))
# #    fitPdf.plotOn(xframe, ROOT.RooFit.LineColor(2))
#     yframe.Draw()

#     can.Print(prefix+'_data%2.1feV_localMinNLL.pdf'%(nuMass1))
#     can.Print(prefix+'_data%2.1feV_localMinNLL.pdf]'%(nuMass1))




    # dataPdf.plotOn(xframe2, ROOT.RooFit.Range('fitRange'),
    #                     ROOT.RooFit.DrawOption("F"),
    #                     ROOT.RooFit.FillColor(ROOT.kOrange),
    #                     ROOT.RooFit.MoveToBack() )

    # fitdata.plotOn(yframe)
    # fitPdf3.plotOn(xframe, ROOT.RooFit.ProjectionRange('fitRange'), 
    #         ROOT.RooFit.Normalization(fitNevts, ROOT.RooAbsReal.NumEvent),
    #         ROOT.RooFit.LineColor(2))
    # fitPdf.plotOn(xframe, ROOT.RooFit.ProjectionRange('fitRange'), 
    #         ROOT.RooFit.Normalization(fitNevts, ROOT.RooAbsReal.NumEvent),
    #         ROOT.RooFit.LineColor(1))
    # dataPdf.plotOn(yframe2, ROOT.RooFit.Range('fitRange'),
    #                     ROOT.RooFit.DrawOption("F"),
    #                     ROOT.RooFit.FillColor(ROOT.kOrange),
    #                     ROOT.RooFit.MoveToBack() )

#     #dataPdf.plotOn(xframe, ROOT.RooFit.Range('fitRange'), ROOT.RooFit.LineColor(4))

#     # can.cd(2)
#     # yframe = nuEnergy.frame(ROOT.RooFit.Title("Test"), 
#     #                         ROOT.RooFit.Bins(100))
#     # yframe.GetYaxis().SetTitleOffset(1.4)
#     # #data.plotOn(yframe)
#     # #deltaT.setVal(bestfitDT)
#     # #fitPdf.plotOn(yframe, ROOT.RooFit.Range('fitRange'))
#     # #deltaT.setVal(0.2)
#     # #fitPdf.plotOn(yframe)
#     # ##deltaT.setVal(0.1)
#     # #fitPdf.plotOn(yframe, ROOT.RooFit.LineColor(2))
#     # yframe.Draw()
#     # data.plotOn(yframe)
#     # # fitPdf3.plotOn(yframe, #ROOT.RooFit.Range('fitRange'), 
#     # #                 ROOT.RooFit.LineColor(2))
#     # # fitPdf.plotOn(yframe, #ROOT.RooFit.Range('fitRange'), 
#     # #                 ROOT.RooFit.LineColor(1))
#     # # yframe.Draw('same')

# #    if dataFile:
# #        dataFile.Close()
# #    if pdfFile:
# #        pdfFile.Close()
# # deltaT.setVal(-0.01)
# # print('nll: ', nll.getVal())

# # fitPdf.fitTo(data, Range('fitTRange'), ROOT.RooFit.PrintEvalErrors(10))
# # nll = ROOT.RooNLLVar("nll", "nll", fitPdf, data,
# #             ROOT.RooFit.Range('fitRange')) #, ROOT.RooFit.Extended())
# # print('nll: ', nll.getVal())
# #nll = fitPdf.createNLL(data)

# # TandE = ROOT.RooArgSet(nuTime, nuEnergy)
# # integr = fitPdf.createIntegral(TandE, ROOT.RooFit.NormSet(TandE), 
# #                 ROOT.RooFit.Range('fitRange'))
# # print('integr: ', integr.getVal())

# #plotPdf = fitPdf

# ###################################
# # 2D Fit to the data set

# # ###################################
# # Draw the PDF
# # can = ROOT.TCanvas("2d fit", "2d fit", 1200, 400)
# # can.Divide(3)
# # can.cd(1)
# # xframe = nuTime.frame(ROOT.RooFit.Title(
# #             "Low statistics histogram pdf"), ROOT.RooFit.Bins(100))
# # xframe.GetYaxis().SetTitleOffset(1.4)
# # #data.plotOnXY(xframe, ROOT.RooFit.YVar(nuEnergy))
# # data.plotOn(xframe)
# # #plotPdf.plotOn(xframe, ROOT.RooFit.Range('fitRange'))
# # xframe.Draw()

# # can.cd(2)
# # yframe = nuEnergy.frame(ROOT.RooFit.Title("Test"), 
# #                         ROOT.RooFit.Bins(100))
# # yframe.GetYaxis().SetTitleOffset(1.4)
# # data.plotOn(yframe)
# # #plotPdf.plotOn(yframe, ROOT.RooFit.Range('fitRange'))
# # yframe.Draw()

# # # can.cd(3)
# # # nllframe = deltaT.frame(ROOT.RooFit.Title("deltaT profile"), 
# # #                         ROOT.RooFit.Bins(10))
# # # nll.plotOn(nllframe, ROOT.RooFit.LineStyle(ROOT.kDashed),
# # #                ROOT.RooFit.LineColor(ROOT.kRed),
# # # #               ROOT.RooFit.PrintEvalErrors(-1),
# # #                ROOT.RooFit.ShiftToZero(),
# # #                ROOT.RooFit.EvalErrorValue(nll.getVal()+10) )
# # # nllframe.SetMaximum(15.)
# # # nllframe.SetMinimum(0.)

# # can.Print('xyfit.pdf[', 'pdf')
# # can.Print('xyfit.pdf')

# # # Draw the data
# # # can.cd(1)
# # # ROOT.gPad.SetLeftMargin(0.15)
# # # inputDataHist.GetZaxis().SetTitleOffset(1.8)
# # # inputDataHist.Draw("surf")

# # # can.cd(2)
# # # ROOT.gPad.SetLeftMargin(0.15)
# # # dataPdf.GetZaxis().SetTitleOffset(1.8)
# # # dataPdf.Draw("surf")

# # # can.cd(3)
# # # fitPdf.GetZaxis().SetTitleOffset(1.8)
# # # fitPdf.Draw("surf")

# # # #frame.Draw()
# # #can.Print("xyfit.pdf")
# # can.Print('xyfit.pdf]', 'pdf')

# # shiftedNuT = ROOT.RooFormulaVar('shiftedNuT', 'shiftedNuT', '@0+@1', 
# #                 ROOT.RooArgList(nuTime, deltaT) )
# # deltaT = ROOT.RooRealVar('deltaT', 'deltaT', 0.0, -0.05, 0.15)
# # data = dataPdf.generate(ROOT.RooArgSet(nuTime, nuEnergy), fitNevts).binnedClone()

# # NevtMin = fitNevts - math.sqrt(fitNevts)*4
# # NevtMax = fitNevts + math.sqrt(fitNevts)*4
# # Nevt = ROOT.RooRealVar('Nevt', 'Nevt', fitNevts, NevtMin, NevtMax)

#     # fitPdf = ROOT.RooHistPdf("fitPdf", "fitPdf", 
#     #             ROOT.RooArgList(shiftedNuT, nuEnergy), 
#     #             ROOT.RooArgList(nuTime, nuEnergy), 
#     #             pdfHist, 1)



#     # nll = ROOT.RooNLLVar("nll", "nll", fitPdf, data, 
#     #         ROOT.RooFit.Range('fitRange')) #ROOT.RooFit.Extended())
#     # nllVal = nll.getVal()

    # nllList, tScan = [], []
    # dTexp = (nuMass1 - nuMass2) * 0.01
    # dTstep = 0.0005
    # NScans = 20
    # locMinNLL, bestfitDT = 100.0, 1.0
    # for iScan in range(NScans):
    #     dT = dTexp - (2*iScan-NScans)*dTstep/2
    #     #dT = -0.005 - iScan*0.0015
    #     deltaT.setVal(dT)

    #     # shift the nuTime in the toy dataset
    #     for iEvt in range(data.sumEntries()):
    #         argSet = data.get(iEvt)
    #         timeTmp = argSet.getRealValue('nuTime')
    #         argSet.setRealValue('nuTime', timeTmp+dT)

    #     nll = ROOT.RooNLLVar("nll", "nll", fitPdf, data, 
    #           ROOT.RooFit.Range('fitRange')) #ROOT.RooFit.Extended())
    #     nllVal = nll.getVal()

    #     if math.isnan(nllVal) != True:
    #         tScan.append( dT )
    #         nllList.append( nllVal )
    #         if locMinNLL>nllVal:
    #             locMinNLL = nllVal
    #             bestfitDT = dT
    
    # print('nll: ', nllList)

#     # file = open(prefix+'_data%2.1feV_pdf%2.1feV.txt'%(nuMass1, nuMass2), 'w')
#     # for k in range(len(nllList)):
#     #     file.write('%8.4f  %10.3f\n'%(tScan[k], nllList[k]))
#     # file.close()

#     # fig2 = plt.figure(figsize=(5,5))
#     # ax = fig2.add_subplot()
#     # ec = plt.plot( tScan, nllList, 'ro')
#     # ax.set_xlabel('$\delta_{T}$ (MeV)')
#     # plt.show()

        # nll = ROOT.RooNLLVar("nll", "nll", fitPdf, fitdata, 
        #                 ROOT.RooFit.Range('nllRange')) #ROOT.RooFit.Extended())
        # nllVal = nll.getVal()

        # #compDataSet(fitdata, reducedData)
        # plotx = np.arange(-0.2, 0.2)
        # plt.plot(oldTList, energyList,  "o", ms=3, label="Sinnock 1969")
        # plt.plot(newTList, energyList,  "o", ms=3, label="Sinnock 1969",color='r')
        # plt.xlabel('$\delta_{T}$ (s)')
        # plt.ylabel('Energy')

    testArgSet = ROOT.RooArgSet(nuTime, nuEnergy)

    testArgSet.setRealValue('nuTime', 0.225)
    testArgSet.setRealValue('nuEnergy', 1.8)
    print('fitPdf2, ', fitPdf.getLogVal( testArgSet ), math.exp(fitPdf.getLogVal( testArgSet )))
    
    integr = fitPdf.createIntegral(testArgSet, ROOT.RooFit.NormSet(testArgSet), 
                ROOT.RooFit.Range('nllRange'))
    print('integr: ', integr.getVal())

    testArgSet.setRealValue('nuTime', 0.013)
    testArgSet.setRealValue('nuEnergy', 2.96)
    print('fitPdf, ', fitPdf.getLogVal( testArgSet ), math.exp(fitPdf.getLogVal( testArgSet )) )