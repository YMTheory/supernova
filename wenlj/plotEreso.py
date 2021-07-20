#!/usr/bin/env python
import ROOT
import math
import sys
import string
from array import array

MaxRotatedEnergy = 12000.

canvas = None
FileNames = { }
ssScintMean = None 
ssIonizMean = None 
msScintMean = None 
msIonizMean = None 

def DivideWithErrors(x, y):
    # If x and y have uncorrelated errors, return their quotient (with errors).
    return (x[0]/y[0],
            math.sqrt(pow(x[1]/y[0], 2) +
                      pow(x[0]*y[1] / pow(y[0], 2), 2)))

def getValue(strline):
    values = strline.split()
    mean = string.atof(values[0])
    err = string.atof(values[1])
    return (mean, err)

def calibE(E, par):
    return (par[0] + E*par[1] + E*E*par[2])

def DoFit(dataSet, PeakPosGuess, SigmaGuess, ErfcFracGuess, ID = None):
    # Using input guesses, pick a fit window and refine those guesses.
    # If ID is given, use it to save a plot.
    RotatedEnergy = dataSet.get()['RotatedEnergy']

    # Note, the bounds used here imply resolution can never be worse than 30% at the peak.
    mean = ROOT.RooRealVar('Mean', 'Mean', PeakPosGuess, PeakPosGuess - SigmaGuess, PeakPosGuess + SigmaGuess)
    sigma = ROOT.RooRealVar('Sigma', 'Sigma', SigmaGuess, 0., 0.3*PeakPosGuess)
    erfc_coeff = ROOT.RooRealVar('coeff', 'Erfc coeff', ErfcFracGuess, 0., 1.)
    GaussianComponent = ROOT.RooGaussian('GaussianComponent', 'GaussianComponent',
                                         RotatedEnergy, mean, sigma)
    ErfcComponent = ROOT.RooGenericPdf('ErfcComponent', 'ErfcComponent',
                                       'TMath::Erfc((@0 - @1) / (TMath::Sqrt(2)*@2))',
                                       ROOT.RooArgList(RotatedEnergy, mean, sigma))
    FitFunction = ROOT.RooAddPdf('FitFunction', 'FitFunction',
                                 GaussianComponent, ErfcComponent, erfc_coeff)
    FitFunction.fitTo(dataSet, ROOT.RooFit.Range(PeakPosGuess - 3*SigmaGuess, PeakPosGuess + 2.5*SigmaGuess))

    if ID:
        # ID[0] should identify the fit as single-site or multi-site.
        # ID[1] should be a formatted identifier of theta.

        # The plot should be wider than the fit window, so we can see what wasn't fit.
        frame = RotatedEnergy.frame(ROOT.RooFit.Bins(100), # ReadME: TO BE ADJUSTED
                                    ROOT.RooFit.Range(PeakPosGuess - 6*SigmaGuess,
                                                      PeakPosGuess + 6*SigmaGuess))
        dataSet.plotOn(frame,
                       ROOT.RooFit.DrawOption('ZP'),
                       ROOT.RooFit.MarkerStyle(20),
                       ROOT.RooFit.MarkerSize(0.5))
        FitFunction.plotOn(frame)
        FitFunction.plotOn(frame,
                           ROOT.RooFit.Components("ErfcComponent"),
                           ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2))
        FitFunction.plotOn(frame,
                           ROOT.RooFit.Components("GaussianComponent"),
                           ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2))
        FitFunction.paramOn(frame,
                            ROOT.RooFit.Layout(0.6))
        frame.SetTitle(ID[0] + ' rotated spectrum (Theta = ' + ID[1] + ' rad)')
        frame.SetXTitle('Rotated Energy (uncalibrated), (Theta = ' + ID[1] + ' rad)')
        frame.Draw()
        canvas.Print(FileNames[ID[0][:2]])

    return ((mean.getVal(), mean.getError()),
            (sigma.getVal(), sigma.getError()),
            (erfc_coeff.getVal(), sigma.getError()))

def GetResolutionForTheta(ID, EnergyList2D, Theta):
    # Build rotated energy dataset for this theta
    RotatedEnergy = ROOT.RooRealVar('RotatedEnergy', 'RotatedEnergy', 0., MaxRotatedEnergy)
    ArgSet = ROOT.RooArgSet(RotatedEnergy)
    dataSet = ROOT.RooDataSet('dataSet', 'dataset', ArgSet)
    for Energy2D in EnergyList2D:
        RotatedEnergy.setVal(Energy2D[0]*math.cos(Theta) + Energy2D[1]*math.sin(Theta))
        dataSet.add(ArgSet)

    # Find the peak.
#    InitPeakPos = FindPeak(ID, dataSet, Theta)
    InitPeakPos = 0
    InitSigma= 0
    global ssScintMean
    global ssIonizMean
    global msScintMean
    global msIonizMean

    if ID=='SS Scint':  # README: initial value is important for a good fit
        InitPeakPos = 2400.
        #InitPeakPos = 6550.
        InitSigma = 0.03 * InitPeakPos
    if ID=='SS Ioniz': 
        InitPeakPos = 3800.
        InitSigma = 0.05 * InitPeakPos
    if ID=='MS Scint': 
        InitPeakPos = 2600.
        #InitPeakPos = 7000.
        InitSigma = 0.065 * InitPeakPos
    if ID=='MS Ioniz': 
        InitPeakPos = 2600.
        InitSigma = 0.04 * InitPeakPos
    if ID=='SS':
        InitPeakPos = ssIonizMean[0]*ROOT.TMath.Cos(Theta) + ssScintMean[0]*ROOT.TMath.Sin(Theta)
        InitSigma = (0.01 + 0.1*ROOT.TMath.Power(Theta-0.9, 2)) * InitPeakPos
    if ID=='MS':
        InitPeakPos = msIonizMean[0]*ROOT.TMath.Cos(Theta) + msScintMean[0]*ROOT.TMath.Sin(Theta)
        InitSigma = (0.018 + 0.1*ROOT.TMath.Power(Theta-0.9, 2)) * InitPeakPos

    # Initialize guesses.  (Uncertainties are bogus.)
    PeakPos = (InitPeakPos, 0.1*InitPeakPos)
    Sigma = (InitSigma, 0.01*InitPeakPos)
    ErfcFrac = (0.1, 0.01)

    # Do a series of intermediate fits, to improve guesses.
    # This turns out to be necessary to get smooth results in some cases.
    for NumFit in range(2):
        PeakPos, Sigma, ErfcFrac = DoFit(dataSet, PeakPos[0], Sigma[0], ErfcFrac[0])

    # Final fit.  Save the results.
    PeakPos, Sigma, ErfcFrac = DoFit(dataSet,
                                     PeakPos[0],
                                     Sigma[0],
                                     ErfcFrac[0],
                                     ID = (ID, '%.5f' % Theta))

    Resolution = DivideWithErrors(Sigma, PeakPos)
    return Resolution, PeakPos

if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setSilentMode(True)
    #ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(1)

    if len(sys.argv) < 3:
        print("Usage: %s [datafilename] [eff] [noise]"  % (sys.argv[0]))
        sys.exit(1)

    # light collection efficiency
    eff = string.atof( sys.argv[2] )
    print 'eff', eff

    noise = string.atof( sys.argv[3] )
    print 'noise', noise

    #dataname = '12mm_rotated_cal'
    dataname = sys.argv[1]
    print 'dataname', dataname 

    chain = ROOT.TChain('data')
    chain.Add('plot/%s.root'%(dataname))
    #for i in range(100):
    #    chain.Add('plot/%s.root'%(dataname, dataname, i))

    Nphoton_ini = array('d', [0.])
    Ncharge_cor = array('d', [0.])

    chain.SetBranchAddress('Nphoton_init', Nphoton_ini)
    chain.SetBranchAddress('Ncharge_cor', Ncharge_cor)

    nEntries = chain.GetEntries()
    print 'nEntries: ', nEntries

    # random generator
    rand = None

    hist = {}
    hist['scintEraw'] = ROOT.TH1F('scintEraw', '', 600, 0, 3500)
    hist['ionizEraw'] = ROOT.TH1F('ionizEraw', '', 600, 0, 3500)
#    hist['scintEtrue'] = ROOT.TH1F('scintEtrue', '', 600, 0, 3500)
#    hist['ionizEtrue'] = ROOT.TH1F('ionizEtrue', '', 500, 0, 1e-3)
    hist['scintE'] = ROOT.TH1F('scintE', '', 600, 0, 3500)
    hist['ionizE'] = ROOT.TH1F('ionizE', '', 600, 0, 3500)
    hist['anti'] = ROOT.TH2F('anti', '', 600, 0, 3500, 600, 0, 3500 )

    hist['Ncharge'] = ROOT.TH1D('Ncharge', '', 500, 0, 100000)
    hist['Nphoton'] = ROOT.TH1D('Nphoton', '', 500, 0, 200000)

    xTitle = { 'scintE': 'scintillation energy (keV)',
               'ionizE': 'ionization energy (keV)',
               'scintEraw': 'scintillation energy (smeared)',
               'ionizEraw': 'ionization energy (smeared)',
               'scintEtrue': 'true scintillation energy (keV)',
               'ionizEtrue': 'ionization thermal electron total energy',
               'Ncharge': 'ionization thermal electron number',
               'Nphoton': 'scintillation photon number',
               'anti':   'ionization energy (keV)'
             }
    yTitle = { 'scintE': 'counts',
               'ionizE': 'counts',
               'scintEraw': 'counts',
               'ionizEraw': 'counts',
               'scintEtrue': 'counts',
               'ionizEtrue': 'counts',
               'Ncharge': 'counts',
               'Nphoton': 'counts',
               'anti':   'scinitllation energy (keV)'

             }
    # set histogram titles
    for name in hist.iterkeys():
        hist[name].GetXaxis().SetTitle( xTitle[name] )
        hist[name].GetYaxis().SetTitle( yTitle[name] )

    # calibration parameters from Eraw
    # for GammaEvents
#    pS = [ -26.068, 0.9258, 5.36521e-06 ]
#    pI = [ 36.5289, 0.9559, -6.74776e-06 ]

    EnergyList2D_ss = []

    for jentry in range(nEntries):
        chain.GetEntry( jentry )
        if jentry==0:
            seed = int( Nphoton_ini[0]/100 )
            print 'seed: ', seed
            rand = ROOT.TRandom3(seed)

        if jentry%1000==0:
            print jentry

        hist['Nphoton'].Fill(Nphoton_ini[0])
        hist['Ncharge'].Fill(Ncharge_cor[0])

        # light collection effect
        # assume 20% light collection efficiency
        #Nphotons = rand.Poisson( eff*Nphoton_ini[0] )

        # energy smearing due to noise?
        #Nnoise = noise*Nphotons
        #Nphotons = Nphotons + rand.Gaus( 0, Nnoise) 

        # raw energy
        scintE = Nphoton_ini[0]/(40.*eff) # FIXME: 40 photons/keV, rough light yield
        ionizE = Ncharge_cor[0]/30. # FIXME: 30 electrons/keV, rough calib.
        hist['scintEraw'].Fill( scintE )
        hist['ionizEraw'].Fill( ionizE )

        # calibrate the energy
        #scintE = calibE(scintE, pS)
        #ionizE = calibE(ionizE, pI)
        hist['scintE'].Fill( scintE )
        hist['ionizE'].Fill( ionizE )
        hist['anti'].Fill( ionizE, scintE )

        EnergyList2D_ss.append( (ionizE, scintE) )

#    ######################################################
    global canvas
    global FileNames
    FileNames['SS'] = 'plot/' + dataname + '_lightEff%5.3f_noise%5.3f.pdf'%(eff, noise)
#    FileNames['MS'] = dataname + '_lightEff%5.3f.pdf'%eff
    canvas = ROOT.TCanvas()
    for name in FileNames.itervalues(): canvas.Print(name + '[')

    # determine the optimal rotation angle and energy resolution
    ThetaToTry = [0.4 + 0.05*i for i in range(20)]
    #ThetaToTry = [0.1 + 0.01*i for i in range(90)]
    graph_ss = ROOT.TGraphErrors()
#    graph_ms = ROOT.TGraphErrors()

    BestTheta_ss = None
    BestRes_ss = None
#    BestTheta_ms = None
#    BestRes_ms = None

    global ssScintMean
    global ssIonizMean
    #global msScintMean
    #global msIonizMean
    ssScintResolution, ssScintMean = GetResolutionForTheta('SS Scint', EnergyList2D_ss, ROOT.TMath.Pi()/2)
    ssIonizResolution, ssIonizMean = GetResolutionForTheta('SS Ioniz', EnergyList2D_ss, 0.)
    #msScintResolution, msScintMean = GetResolutionForTheta('MS Scint', EnergyList2D_ms, ROOT.TMath.Pi()/2)
    #msIonizResolution, msIonizMean = GetResolutionForTheta('MS Ioniz', EnergyList2D_ms, 0.)

    print ssIonizMean, ssScintMean
    print ssIonizMean[0], ssScintMean[0]
    #print msIonizMean, msScintMean
    #print msIonizMean[0], msScintMean[0]

    for TestTheta in ThetaToTry:
#        print "TestTheta ", TestTheta
        ssResolution, _ = GetResolutionForTheta('SS', EnergyList2D_ss, TestTheta)
        #msResolution, _ = GetResolutionForTheta('MS', EnergyList2D_ms, TestTheta)

        graph_ss.Set(graph_ss.GetN()+1)
        graph_ss.SetPoint(graph_ss.GetN()-1, TestTheta, ssResolution[0])
        graph_ss.SetPointError(graph_ss.GetN()-1, 0., ssResolution[1])
#        graph_ms.Set(graph_ms.GetN()+1)
#        graph_ms.SetPoint(graph_ms.GetN()-1, TestTheta, msResolution[0])
#        graph_ms.SetPointError(graph_ms.GetN()-1, 0., msResolution[1])

        if BestRes_ss == None or ssResolution[0] < BestRes_ss:
            BestRes_ss = ssResolution[0]
            BestTheta_ss = TestTheta
#        if BestRes_ms == None or msResolution[0] < BestRes_ms:
#            BestRes_ms = msResolution[0]
#            BestTheta_ms = TestTheta
#
    # ToDO:  Consider trying to pick a fit window more intelligently.
    ssPolynomialFunc = ROOT.TF1("ssFit", "[0] + [1]*(x-[2])**2", BestTheta_ss - 0.2, BestTheta_ss + 0.2)
    ssPolynomialFunc.SetParameters(BestRes_ss, 4.0, BestTheta_ss)
    graph_ss.Fit(ssPolynomialFunc, "ER")
#    msPolynomialFunc = ROOT.TF1("msFit", "[0] + [1]*(x-[2])**2", BestTheta_ms - 0.07, BestTheta_ms + 0.07)
#    msPolynomialFunc.SetParameters(BestRes_ms, 4.0, BestTheta_ms)
#    graph_ms.Fit(msPolynomialFunc, "ER")

    ROOT.gStyle.SetOptFit(11)
    graph_ss.SetTitle('SS Resolution vs Theta')
    graph_ss.GetXaxis().SetTitle('Theta (radians)')
    graph_ss.GetYaxis().SetTitle('Resolution (at thorium peak)')
    graph_ss.Draw("A*")
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs(FileNames['SS'])
#    graph_ms.SetTitle('MS Resolution vs Theta')
#    graph_ms.GetXaxis().SetTitle('Theta (radians)')
#    graph_ms.GetYaxis().SetTitle('Resolution (at thorium peak)')
#    graph_ms.Draw("A*")
#    canvas.Modified()
#    canvas.Update()
#    canvas.SaveAs(FileNames['MS'])
#
    ssResolution, ssMean = GetResolutionForTheta('SS', EnergyList2D_ss, ssPolynomialFunc.GetParameter(2))
##    msResolution, msMean = GetResolutionForTheta('MS', EnergyList2D_ms, msPolynomialFunc.GetParameter(2))
#
    for name in FileNames.itervalues(): canvas.Print(name + ']')

    result = { 'Theta_ss' : (ssPolynomialFunc.GetParameter(2), ssPolynomialFunc.GetParError(2)),
             'PeakPos_ss' : ssMean,
             'Resolution_ss' : ssResolution,
#             'Theta_ms' : (msPolynomialFunc.GetParameter(2), msPolynomialFunc.GetParError(2)),
#             'PeakPos_ms' : msMean,
#             'Resolution_ms' : msResolution,
             'ssScintPeak' : ssScintMean,
             'ssScintRes' : ssScintResolution,
             'ssIonizPeak' : ssIonizMean,
             'ssIonizRes' : ssIonizResolution,
#             'msScintPeak' : msScintMean,
#             'msScintRes' : msScintResolution,
#             'msIonizPeak' : msIonizMean,
#             'msIonizRes' : msIonizResolution }
           }
#
#    #######################################################
#    c1 = ROOT.TCanvas('c1', 'c1', 1200, 500)
#    c1.Divide(2,1)
#    c1.cd(1)
#    hist['scintEtrue'].Draw()
#    hist['scintEtrue'].Fit('gaus', '+')
#    c1.cd(2)
#    hist['ionizEtrue'].Draw()
#    hist['ionizEtrue'].Fit('gaus', '+')
#    c1.Print('gamma/%s_Etrue.png'%dataname)
#
#    c2 = ROOT.TCanvas('c2', 'c2', 1200, 500)
#    c2.Divide(2,1)
#    c2.cd(1)
#    hist['scintEraw'].Draw()
#    hist['scintEraw'].Fit('gaus', '+')
#    c2.cd(2)
#    hist['ionizEraw'].Draw()
#    hist['ionizEraw'].Fit('gaus', '+')
#    c2.Print('gamma/%s_Eraw_lightEff%5.3f_noise%5.3f.png'%(dataname, eff, noise))
#
#    c5 = ROOT.TCanvas('c5', 'c5', 1200, 500)
#    c5.Divide(2,1)
#    c5.cd(1)
#    hist['Nphoton_ini'].Draw()
#    hist['Nphoton_ini'].Fit('gaus', '+')
#    c5.cd(2)
#    hist['TECounter'].Draw()
#    hist['TECounter'].Fit('gaus', '+')
#    c5.Print('gamma/%s_SignalRaw.png'%dataname)
#
#    c3 = ROOT.TCanvas('c3', 'c3', 600, 600)
#    hist['anti'].Draw('COLZ')
#    c3.Print('gamma/%s_Anti_lightEff%5.3f_noise%5.3f.png'%(dataname,eff,noise))
#
#    c4 = ROOT.TCanvas('c4', 'c4', 600, 500)
#    hist['ionizE'].Draw()
#    hist['scintE'].Draw('same')
#    hist['scintE'].SetLineColor(2)
#    c4.Print('gamma/%s_calibE_lightEff%5.3f_noise%5.3f.png'%(dataname, eff, noise))
#
    output = ROOT.TFile('plot/'+dataname+'_lightEff%5.3f_noise%5.3f.root'%(eff, noise), 'recreate')
    for name in hist.iterkeys():
        hist[name].Write()
    output.Close()
#
    file = open('plot/'+dataname+'_lightEff%5.3f_noise%5.3f.txt'%(eff, noise), 'w')
    resNames = [ 'Theta_ss', 'PeakPos_ss', 'Resolution_ss',
                 'ssScintPeak', 'ssScintRes', 'ssIonizPeak', 'ssIonizRes' ]
    for name in resNames:
        value = result[name]
        if type(value) == type(1.0):
            file.write( str(float(value))+"\n")
        if type(value) == type(()):
            file.write( str(float(value[0]))+"  "+str(float(value[1]))+"\n"  )
    file.close()

#    raw_input()
