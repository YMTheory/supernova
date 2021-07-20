#!/usr/bin python
#wrap_python.sh plotGamma.py EFFICIENCY NOISE
import ROOT
import math
import sys
import string
import os
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

def calibScintE(E, par):
    return (par[0] + E*par[1] + E*E*par[2])

def calibIonizE(E, par):
    return (par[0] + E*par[1] + E*E*par[2])
    #return E/(ROOT.TMath.Exp(-par[1]/par[0]))

def DoFit(dataHist, PeakPosGuess, SigmaGuess, ErfcFracGuess, ID = None):
    # Using input guesses, pick a fit window and refine those guesses.
    # If ID is given, use it to save a plot.
    RotatedEnergy = dataHist.get()['RotatedEnergy']

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
    FitFunction.fitTo(dataHist, ROOT.RooFit.Range(PeakPosGuess - 5*SigmaGuess, PeakPosGuess + 3.5*SigmaGuess))
    #FitFunction.fitTo(dataHist, ROOT.RooFit.Range(PeakPosGuess - 1.0*SigmaGuess, PeakPosGuess + 3.5*SigmaGuess))

    if ID:
        # ID[0] should identify the fit as single-site or multi-site.
        # ID[1] should be a formatted identifier of theta.

        # The plot should be wider than the fit window, so we can see what wasn't fit.
        frame = RotatedEnergy.frame(ROOT.RooFit.Bins(200),
                                    ROOT.RooFit.Range(PeakPosGuess - 6*SigmaGuess,
                                                      PeakPosGuess + 6*SigmaGuess))
        dataHist.plotOn(frame,
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
        frame.SetXTitle('Rotated Energy (uncalibrated)')
        frame.Draw()
        canvas.Print(FileNames[ID[0][:2]])

    return ((mean.getVal(), mean.getError()),
            (sigma.getVal(), sigma.getError()),
            (erfc_coeff.getVal(), sigma.getError()))

def GetResolutionForTheta(ID, EnergyList2D, Theta):
    # Build rotated energy dataset for this theta
    RotatedEnergy = ROOT.RooRealVar('RotatedEnergy', 'RotatedEnergy', 0., MaxRotatedEnergy)
#    ArgSet = ROOT.RooArgSet(RotatedEnergy)
#    dataSet = ROOT.RooDataSet('dataSet', 'dataset', ArgSet)
    hist = ROOT.TH1D('hist', '', 1000, 0, 5000)
    for Energy2D in EnergyList2D:
        hist.Fill( Energy2D[0]*math.cos(Theta) + Energy2D[1]*math.sin(Theta) )
#        RotatedEnergy.setVal(Energy2D[0]*math.cos(Theta) + Energy2D[1]*math.sin(Theta))
#        dataSet.add(ArgSet)

    # RooDataHist
    dataHist = ROOT.RooDataHist('dataSet', 'dataSet',
               ROOT.RooArgList(RotatedEnergy), hist, 1.0)

    # Find the peak.
#    InitPeakPos = FindPeak(ID, dataSet, Theta)
    InitPeakPos = 0
    InitSigma= 0
    global ssScintMean
    global ssIonizMean
    global msScintMean
    global msIonizMean

    if ID=='SS Scint': 
        InitPeakPos = hist.GetMean()
        #InitPeakPos = hist.GetBinCenter(hist.GetMaximumBin())
        #InitPeakPos = 6550.
        InitSigma = hist.GetRMS()
    if ID=='SS Ioniz': 
        InitPeakPos = hist.GetMean()
        #InitPeakPos = hist.GetBinCenter(hist.GetMaximumBin())
        InitSigma = hist.GetRMS()
    if ID=='MS Scint': 
        InitPeakPos = 2500.
        #InitPeakPos = 7000.
        InitSigma = 0.065 * InitPeakPos
    if ID=='MS Ioniz': 
        InitPeakPos = 2500.
        InitSigma = 0.04 * InitPeakPos
    if ID=='SS':
        InitPeakPos = ssIonizMean[0]*ROOT.TMath.Cos(Theta) + ssScintMean[0]*ROOT.TMath.Sin(Theta)
        InitSigma = (0.01 + 0.1*ROOT.TMath.Power(Theta-0.9, 2)) * InitPeakPos
    if ID=='MS':
        InitPeakPos = msIonizMean[0]*ROOT.TMath.Cos(Theta) + msScintMean[0]*ROOT.TMath.Sin(Theta)
        InitSigma = (0.01 + 0.1*ROOT.TMath.Power(Theta-0.9, 2)) * InitPeakPos

    hist.Delete()

    # Initialize guesses.  (Uncertainties are bogus.)
    PeakPos = (InitPeakPos, 0.1*InitPeakPos)
    Sigma = (InitSigma, 0.01*InitPeakPos)
    ErfcFrac = (0.3, 0.01)

    # Do a series of intermediate fits, to improve guesses.
    # This turns out to be necessary to get smooth results in some cases.
    for NumFit in range(2):
        PeakPos, Sigma, ErfcFrac = DoFit(dataHist, PeakPos[0], Sigma[0], ErfcFrac[0])

    # Final fit.  Save the results.
    PeakPos, Sigma, ErfcFrac = DoFit(dataHist,
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
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetPalette(1)
    ROOT.gInterpreter.ExecuteMacro("/junofs/users/wenlj/nEXO/resolution/style.cc");

    if len(sys.argv) < 3:
        print("Usage: %s [type] [runlist]"  % (sys.argv[0]))
        print("[type] should be either 'list' or 'run'")
        sys.exit(1)

    # light collection efficiency
    eff = string.atof( sys.argv[1] )
    print 'eff', eff

    ldarkrate = string.atof( sys.argv[2] )
    print 'dark rate in light channel', ldarkrate

    lnoise = string.atof( sys.argv[3] )
    print 'noise of light channel', lnoise

    lam = string.atof( sys.argv[4] )
    print 'lambda of light channel', lam

    enoise = string.atof( sys.argv[5] )
    print 'noise of charge channel', enoise

    diff = string.atof( sys.argv[6] )
    print 'diffusion', diff 

    pitch = string.atof( sys.argv[7] )
    print 'pitch size', pitch 

    energy = string.atof( sys.argv[8] )
    print 'energy of electron', energy

    efield = string.atof( sys.argv[9] )
    print 'electronic field', efield

    seed = string.atoi( sys.argv[10] )
    print 'random seed', seed

    dataDir = '/junofs/users/wenlj/nEXO/resolution/'
    dataname = 'bb0n_%.0f_EF%s_eff%s_ldarkrate%s_lnoise%s_lambda%s_enoise%s'%(energy, efield, eff, ldarkrate, lnoise, lam, enoise)

    filename = 'Selected_bb0n_waveformcut_5_2M_diffu_%s_noise_%s_wholevolume_PAD_%smm_%sVcm_event10000_el_result'%(sys.argv[6], sys.argv[5], sys.argv[7], sys.argv[8])
    file = open('%s/%s.txt'%(dataDir, filename))
    lines = file.readlines()
    (meantau, meantauErr) = getValue( lines[0] )
    file.close()
    print 'Check EL: ', meantau, meantauErr

    chain = ROOT.TChain('TreeS')
    for i in range(1):
        filename = '/junofs/users/jiay/nEXO/TMVA/output/Selected_bb0n_waveformcut_5_2M_diffu_%s_noise_%s_wholevolume_PAD_%smm_%sVcm_event10000.root'%(sys.argv[6], sys.argv[5], sys.argv[7], sys.argv[8])
        chain.Add(filename)
    
    GenX = array('d', [0.])
    GenY = array('d', [0.])
    GenZ = array('d', [0.])
    Multiplicity = array('i', [0])
    NumScintGammas = array('d', [0])
    ScintGammaskeV = array('d', 1000000*[0.])
    OPTime = array('d', 1000000*[0.])
    SiPMID = array('i', 1000000*[0])
    TECounter = array('d', [0])
    TEKineticEnergy = array('d', 1000000*[0.])
    TEX = array('d', 1000000*[0.])
    TEY = array('d', 1000000*[0.])
    TEZ = array('d', 1000000*[0.])
    driftTime = array('d', [0])

    #chain.SetBranchAddress('GenX', GenX)
    #chain.SetBranchAddress('GenY', GenY)
    #chain.SetBranchAddress('GenZ', GenZ)
    #chain.SetBranchAddress('multiplicity', Multiplicity)
    chain.SetBranchAddress('NumScint', NumScintGammas)
    #chain.SetBranchAddress('OPEnergy', ScintGammaskeV)
    #chain.SetBranchAddress('OPTime', OPTime)
    #chain.SetBranchAddress('SiPMID', SiPMID)
    chain.SetBranchAddress('NumElectron', TECounter )
    chain.SetBranchAddress('DriftTime', driftTime)
    #chain.SetBranchAddress('TEEnergy', TEKineticEnergy )
    #chain.SetBranchAddress('TEX', TEX )
    #chain.SetBranchAddress('TEY', TEY )
    #chain.SetBranchAddress('TEZ', TEZ )

    nEntries = chain.GetEntries()
    print 'nEntries: ', nEntries

    # create a root file to store tree
    evtTree = ROOT.TTree("evtTree","ibdTree to save truth info.")
    scintE_tree = array('d',[0])
    ionizE_tree = array('d',[0])

    evtTree.Branch('scintE', scintE_tree,'scintE/D')
    evtTree.Branch('ionizE', ionizE_tree,'ionizE/D')

    # random generator
    rand = None

    hist = {}
    hist['prob_ct'] = ROOT.TH1F('prob_ct', '', 10, 1, 11)
    hist['scintEraw'] = ROOT.TH1F('scintEraw', '', 800, 0, 4000)
    hist['ionizEraw'] = ROOT.TH1F('ionizEraw', '', 800, 0, 4000)
    hist['scintEtrue'] = ROOT.TH1F('scintEtrue', '', 500, 0, 200000)
    hist['ionizEtrue'] = ROOT.TH1F('ionizEtrue', '', 500, 0, 1e-3)
    hist['scintE'] = ROOT.TH1F('scintE', '', 800, 0, 4000)
    hist['ionizE'] = ROOT.TH1F('ionizE', '', 800, 0, 4000)
    hist['anti'] = ROOT.TH2F('anti', '', 800, 0, 4000, 800, 0, 4000 )
    hist['purity'] = ROOT.TH2F('purity', '', 200, 0, 1000, 800, 0, 4000 )

    hist['TECounter'] = ROOT.TH1F('TECounter', '', 1000, 0, 200000)
    hist['NumScintGammas'] = ROOT.TH1F('NumScintGammas', '', 500, 0, 200000)

    xTitle = { 'prob_ct': 'number of pe',
               'scintE': 'scintillation energy (keV)',
               'ionizE': 'ionization energy (keV)',
               'scintEraw': 'scintillation energy (smeared)',
               'ionizEraw': 'ionization energy (smeared)',
               'scintEtrue': 'true scintillation energy (keV)',
               'ionizEtrue': 'ionization thermal electron total energy',
               'TECounter': 'ionization thermal electron number',
               'NumScintGammas': 'scintillation photon number',
               'anti':   'ionization energy (keV)',
               'purity':   'drift time (#mus)'
             }
    yTitle = { 'prob_ct': 'prob',
               'scintE': 'counts',
               'ionizE': 'counts',
               'scintEraw': 'counts',
               'ionizEraw': 'counts',
               'scintEtrue': 'counts',
               'ionizEtrue': 'counts',
               'TECounter': 'counts',
               'NumScintGammas': 'counts',
               'anti':   'scinitllation energy (keV)',
               'purity':   'raw ionization energy (keV)'
             }
    # set histogram titles
    for name in hist.iterkeys():
        hist[name].GetXaxis().SetTitle( xTitle[name] )
        hist[name].GetYaxis().SetTitle( yTitle[name] )

    # calibration parameters from Eraw
    # for GammaEvents
#    pS = [ -26.068, 0.9258, 5.36521e-06 ]
#    pI = [ 36.5289, 0.9559, -6.74776e-06 ]
    pS = [ 0., 1., 0. ]
    pI = [ 0., 1., 0. ]

    # define probability function of crosstalk.
    for ibin in range(10):
        prob = ROOT.TMath.Power(lam*(ibin+1),ibin)*ROOT.TMath.Exp(-lam*(ibin+1))/ROOT.TMath.Factorial(ibin+1)
        hist['prob_ct'].SetBinContent(ibin+1,prob)

    EnergyList2D_ss = []

    rand = ROOT.TRandom3(seed)

#    for jentry in range(nEntries):
#        chain.GetEntry( jentry )
#        if Multiplicity[0] > 1:
#           continue
#        if TECounter[0] > 1e6:
#           continue
#        ThermalElectrons = TECounter[0]
#        NfinalThermalElectrons = ThermalElectrons * ROOT.TMath.Exp( driftTime[0]/meantau )
#
#        #Zposition = rand.Rndm()*1259.0
#        #deltaT = Zposition/velocity
#        #tau = rand.Gaus(meantau, sigmatau)
#        #ThermalElectrons - rand.Poisson(ThermalElectrons*(1-ROOT.TMath.Exp(-deltaT/tau)))
#        #NfinalThermalElectrons = ThermalElectrons - rand.Poisson(ThermalElectrons*(1-ROOT.TMath.Exp(-deltaT/tau)))
#        #NfinalThermalElectrons = NfinalThermalElectrons + rand.Gaus(0, enoise)
#
#        ionizE = NfinalThermalElectrons/(1.14831e+05/2500) # FIXME
#        #hist['purity'].Fill( deltaT, ionizE )
#
#    #prof = hist['purity'].ProfileX()
#    #prof.Fit('expo', '+')
#    #fun = prof.GetFunction('expo')
#    #fittedTau = fun.GetParameter(1)
#    #print -1/fittedTau
 
    for jentry in range(nEntries):
        chain.GetEntry( jentry )

        #if Multiplicity[0] > 1:
        #   continue

        if TECounter[0] > 5e5:
           continue

        if jentry%1000==0:
            print jentry

        #if Reflectivity[0]!=ref or AbsLength[0]!=abs:
        #   continue

        hist['NumScintGammas'].Fill(NumScintGammas[0])
        hist['TECounter'].Fill(TECounter[0])

        #for iOP in range(ScintGammaskeV[0]):
        #    hist['scintEtrue'].Fill( ScintGammaskeV[iOP] )

        #for iTE in range(TECounter[0]):
        #    hist['ionizEtrue'].Fill( TEKineticEnergy[iTE] )

        # light collection effect
        Nphotons = rand.Poisson(eff*NumScintGammas[0])

        # add dark noise
        Nphotons = Nphotons + rand.Poisson(ldarkrate)

        NfinalPhotons = 0
        # add correlated noise
        for iphoton in range(Nphotons):
            NfinalPhotons = NfinalPhotons + int(hist['prob_ct'].GetRandom())

        # add electronics noise
        NfinalPhotons = NfinalPhotons + rand.Gaus(0, lnoise)

        ThermalElectrons = TECounter[0]
        ThermalElectrons = ThermalElectrons * ROOT.TMath.Exp( driftTime[0]/meantau )

        # generated positions
        #Zposition = rand.Rndm()*1259.0

        #deltaT = Zposition/velocity
        #tau = rand.Gaus(meantau, sigmatau)

        #ThermalElectrons = TECounter[0] - rand.Poisson(TECounter[0]*(1-ROOT.TMath.Exp(-deltaT/tau)))

        # assume 2keV noise in ionization channel (800 e RMS --> 16keV, from Liang)
        # elecnoise = 800.
        #ThermalElectrons = ThermalElectrons + rand.Gaus(0, enoise)

        hist['scintEtrue'].Fill( NfinalPhotons )

        # raw energy
        scintEscale = 0.0
        scintE = 0.0
        ionizE = 0.0
      
        #################################################################
        ####   diffusion coefficient = 100 cm^2/s                     ###
        ################################################################# 
        if efield == 65:
#           if energy == 100:
#              scintEscale = (1/(1-lam))*(5.50061e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(553.3/100) # FIXME
#
#           if energy == 500:
#              scintEscale = (1/(1-lam))*(2.30960e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.34923e+04/500) # FIXME
#
#           if energy == 1000:
#              scintEscale = (1/(1-lam))*(4.41778e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(3.10762e+04/1000) # FIXME
#
#           if energy == 1500:
#              scintEscale = (1/(1-lam))*(6.51187e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(4.87064e+04/1500) # FIXME
#
#           if energy == 2000:
#              scintEscale = (1/(1-lam))*(8.59939e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(6.63179e+04/2000) # FIXME

           if energy == 2458:
              #noise = 200, threshold = 1000.
              scintEscale = (1/(1-lam))*(1.05165e+05*eff + ldarkrate)/2458
              scintE = NfinalPhotons/scintEscale # FIXME
              #ionizE = ThermalElectrons/(8.25217e+04/2458) # FIXME
              #noise = 50, threshold = 150.
              #ionizE = ThermalElectrons/(8.30392e+04/2458) # FIXME
              #noise = 50, threshold = 250.
              #ionizE = ThermalElectrons/(8.29648e+04/2458) # FIXME
              #noise = 100, threshold = 300.
              #ionizE = ThermalElectrons/(8.29436e+04/2458) # FIXME
              #noise = 100, threshold = 500.
              #ionizE = ThermalElectrons/(8.27822e+04/2458) # FIXME
              #noise = 200, threshold = 600.
              #ionizE = ThermalElectrons/(8.27408e+04/2458) # FIXME
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(8.25073e+04/2458) # FIXME
              #noise = 400, threshold = 1200.
              ionizE = ThermalElectrons/(8.22774e+04/2458) # FIXME

        if efield == 130:
#           if energy == 100:
#              scintEscale = (1/(1-lam))*(4.92790e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1160.0/100) # FIXME
#
#           if energy == 500:
#              scintEscale = (1/(1-lam))*(2.05941e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.63310e+04/500) # FIXME
#
#           if energy == 1000:
#              scintEscale = (1/(1-lam))*(3.94002e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(3.60803e+04/1000) # FIXME
#
#           if energy == 1500:
#              scintEscale = (1/(1-lam))*(5.80589e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(5.59311e+04/1500) # FIXME
#
#           if energy == 2000:
#              scintEscale = (1/(1-lam))*(7.66532e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(7.57492e+04/2000) # FIXME

           if energy == 2458:
              #noise = 200, threshold = 1000.
              scintEscale = (1/(1-lam))*(9.36506e+04*eff + ldarkrate)/2458
              scintE = NfinalPhotons/scintEscale # FIXME
              #ionizE = ThermalElectrons/(9.38896e+04/2458) # FIXME
              #noise = 50, threshold = 150.
              #ionizE = ThermalElectrons/(9.45456e+04/2458) # FIXME
              #noise = 50, threshold = 250.
              #ionizE = ThermalElectrons/(9.44703e+04/2458) # FIXME
              #noise = 100, threshold = 300.
              #ionizE = ThermalElectrons/(9.44491e+04/2458) # FIXME
              #noise = 100, threshold = 500.
              #ionizE = ThermalElectrons/(9.42837e+04/2458) # FIXME
              #noise = 200, threshold = 600.
              #ionizE = ThermalElectrons/(9.42851e+04/2458) # FIXME
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(9.40466e+04/2458) # FIXME
              #noise = 400, threshold = 1200.
              ionizE = ThermalElectrons/(9.41987e+04/2458) # FIXME

        if efield == 200:
#           if energy == 100:
#              scintEscale = (1/(1-lam))*(4.53828e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1618.0/100) # FIXME
#
#           if energy == 500:
#              scintEscale = (1/(1-lam))*(1.89560e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.80809e+04/500) # FIXME
#
#           if energy == 1000:
#              scintEscale = (1/(1-lam))*(3.63586e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(3.92186e+04/1000) # FIXME
#
#           if energy == 1500:
#              scintEscale = (1/(1-lam))*(5.35564e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(6.04923e+04/1500) # FIXME
#
#           if energy == 2000:
#              scintEscale = (1/(1-lam))*(7.06559e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(8.18062e+04/2000) # FIXME

           if energy == 2458:
              #noise = 200, threshold = 1000.
              #scintEscale = (1/(1-lam))*(8.63334e+04*eff + ldarkrate)/2458
              scintEscale = (1/(1-lam))*(5.55e+04*eff + ldarkrate)/2458
              scintE = NfinalPhotons/scintEscale # FIXME
              #ionizE = ThermalElectrons/(1.01257e+05/2458) # FIXME
              #noise = 50, threshold = 150.
              #ionizE = ThermalElectrons/(1.01909e+05/2458) # FIXME
              #noise = 50, threshold = 250.
              #ionizE = ThermalElectrons/(1.01845e+05/2458) # FIXME
              #noise = 100, threshold = 300.
              #ionizE = ThermalElectrons/(1.01802e+05/2458) # FIXME
              #noise = 100, threshold = 500.
              #ionizE = ThermalElectrons/(1.01648e+05/2458) # FIXME
              #noise = 200, threshold = 600.
              #ionizE = ThermalElectrons/(1.01605e+05/2458) # FIXME
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(1.01431e+05/2458) # FIXME
              #noise = 400, threshold = 1200.
              #ionizE = ThermalElectrons/(1.01549e+05/2458) # FIXME
              ionizE = ThermalElectrons/(1.21549e+05/2458) # FIXME

        if efield == 300:
#           if energy == 100:
#              scintEscale = (1/(1-lam))*(4.16092e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(2069.0/100) # FIXME
#
#           if energy == 500:
#              scintEscale = (1/(1-lam))*(1.74798e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.96242e+04/500) # FIXME
#
#           if energy == 1000:
#              scintEscale = (1/(1-lam))*(3.35502e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(4.20717e+04/1000) # FIXME
#        
#           if energy == 1500:
#              scintEscale = (1/(1-lam))*(4.94373e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(6.46345e+04/1500) # FIXME
#
#           if energy == 2000:
#              scintEscale = (1/(1-lam))*(6.52369e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(8.72463e+04/2000) # FIXME

           if energy == 2458:
              #noise = 200, threshold = 1000.
              scintEscale = (1/(1-lam))*(7.96883e+04*eff + ldarkrate)/2458
              scintE = NfinalPhotons/scintEscale # FIXME
              #ionizE = ThermalElectrons/(1.07909e+05/2458) # FIXME
              #noise = 50, threshold = 150.
              #ionizE = ThermalElectrons/(1.08524e+05/2458) # FIXME
              #noise = 50, threshold = 250.
              #ionizE = ThermalElectrons/(1.08464e+05/2458) # FIXME
              #noise = 100, threshold = 300.
              #ionizE = ThermalElectrons/(1.08461e+05/2458) # FIXME
              #noise = 100, threshold = 500.
              #ionizE = ThermalElectrons/(1.08296e+05/2458) # FIXME
              #noise = 200, threshold = 600.
              #ionizE = ThermalElectrons/(1.08249e+05/2458) # FIXME
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(1.08080e+05/2458) # FIXME
              #noise = 400, threshold = 1200.
              ionizE = ThermalElectrons/(1.08181e+05/2458) # FIXME

        if efield == 400:
#           if energy == 100:
#              #noise = 200, threshold = 1000.
#              #scintEscale = (1/(1-lam))*(3.89933e+03*eff + ldarkrate)/100
#              #scintE = NfinalPhotons/scintEscale # FIXME
#              #ionizE = ThermalElectrons/(2385.0/100) # FIXME
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(3.90042e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(3.05195e+03/100)
#           if energy == 200:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(7.13034e+03*eff + ldarkrate)/200
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(7.49478e+03/200)
#           if energy == 300:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(1.02762e+04*eff + ldarkrate)/300
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.20382e+04/300)
#           if energy == 400:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(1.33777e+04*eff + ldarkrate)/400
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.66230e+04/400)
#
#           if energy == 500:
#              #noise = 200, threshold = 1000.
#              #scintEscale = (1/(1-lam))*(1.64747e+04*eff + ldarkrate)/500
#              #scintE = NfinalPhotons/scintEscale # FIXME
#              #ionizE = ThermalElectrons/(2.06856e+04/500) # FIXME
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(1.64588e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(2.13035e+04/500) # FIXME
#
#           if energy == 1000:
#              #noise = 200, threshold = 1000.
#              #scintEscale = (1/(1-lam))*(3.16811e+04*eff + ldarkrate)/1000
#              #scintE = NfinalPhotons/scintEscale # FIXME
#              #ionizE = ThermalElectrons/(4.39943e+04/1000) # FIXME
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(3.15869e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(4.47121e+04/1000) # FIXME
#
#           if energy == 1500:
#              #noise = 200, threshold = 1000.
#              #scintEscale = (1/(1-lam))*(4.67089e+04*eff + ldarkrate)/1500
#              #scintE = NfinalPhotons/scintEscale # FIXME
#              #ionizE = ThermalElectrons/(6.73998e+04/1500) # FIXME
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(4.64815e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(6.83511e+04/1500) # FIXME
#
#           if energy == 2000:
#              #noise = 200, threshold = 1000.
#              #scintEscale = (1/(1-lam))*(6.15934e+04*eff + ldarkrate)/2000
#              #scintE = NfinalPhotons/scintEscale # FIXME
#              #ionizE = ThermalElectrons/(9.09177e+04/2000) # FIXME
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(6.11726e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(9.21894e+04/2000) # FIXME

           if energy == 2458:
              scintEscale = (1/(1-lam))*(7.52264e+04*eff + ldarkrate)/2458
              #scintEscale = (1/(1-lam))*(7.45956e+04*eff + ldarkrate)/2458 #SS
              scintE = NfinalPhotons/scintEscale # FIXME
              #ionizE = ThermalElectrons/(1.12412e+05/2458) # FIXME
              #noise = 50, threshold = 150.
              #ionizE = ThermalElectrons/(1.13058e+05/2458) # FIXME
              #noise = 50, threshold = 250.
              #ionizE = ThermalElectrons/(1.12988e+05/2458) # FIXME
              #noise = 100, threshold = 300.
              #ionizE = ThermalElectrons/(1.12978e+05/2458) # FIXME
              #noise = 100, threshold = 500.
              #ionizE = ThermalElectrons/(1.12825e+05/2458) # FIXME
              #noise = 200, threshold = 600.
              #ionizE = ThermalElectrons/(1.12811e+05/2458) # SS+MS
              #ionizE = ThermalElectrons/(1.14030e+05/2458) # SS
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(1.12576e+05/2458) # FIXME
              #noise = 400, threshold = 1200.
              ionizE = ThermalElectrons/(1.12709e+05/2458) # FIXME
              #noise = 800, threshold = 2400.
              #ionizE = ThermalElectrons/(1.11725e+05/2458) # FIXME
              #ionizE = ThermalElectrons/(1.12784e+05/2458) # SS

        if efield == 500:
#           if energy == 100:
#              scintEscale = (1/(1-lam))*(3.70931e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(2598.0/100) # FIXME
#
#           if energy == 500:
#              scintEscale = (1/(1-lam))*(1.57250e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(2.14687e+04/500) # FIXME
#
#           if energy == 1000:
#              scintEscale = (1/(1-lam))*(3.02724e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(4.54099e+04/1000) # FIXME
#
#           if energy == 1500:
#              scintEscale = (1/(1-lam))*(4.46746e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(6.94577e+04/1500) # FIXME
#
#           if energy == 2000:
#              scintEscale = (1/(1-lam))*(5.89021e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(9.36175e+04/2000) # FIXME

           if energy == 2458:
              scintEscale = (1/(1-lam))*(7.19294e+04*eff + ldarkrate)/2458
              scintE = NfinalPhotons/scintEscale # FIXME
              #ionizE = ThermalElectrons/(1.15727e+05/2458) # FIXME
              #noise = 50, threshold = 150.
              #ionizE = ThermalElectrons/(1.16308e+05/2458) # FIXME
              #noise = 50, threshold = 250.
              #ionizE = ThermalElectrons/(1.16250e+05/2458) # FIXME
              #noise = 100, threshold = 300.
              #ionizE = ThermalElectrons/(1.16230e+05/2458) # FIXME
              #noise = 100, threshold = 500.
              #ionizE = ThermalElectrons/(1.16069e+05/2458) # FIXME
              #noise = 200, threshold = 600.
              #ionizE = ThermalElectrons/(1.16051e+05/2458) # FIXME
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(1.15879e+05/2458) # FIXME
              #noise = 400, threshold = 1200.
              ionizE = ThermalElectrons/(1.15994e+05/2458) # FIXME

        if efield == 750:
#           if energy == 100:
#              scintEscale = (1/(1-lam))*(3.37444e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(2982.0/100) # FIXME
#
#           if energy == 500:
#              scintEscale = (1/(1-lam))*(1.44603e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(2.27699e+04/500) # FIXME
#
#           if energy == 1000:
#              scintEscale = (1/(1-lam))*(2.79183e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(4.78060e+04/1000) # FIXME
#
#           if energy == 1500:
#              scintEscale = (1/(1-lam))*(4.11863e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(7.29848e+04/1500) # FIXME
#
#           if energy == 2000:
#              scintEscale = (1/(1-lam))*(5.43534e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(9.82095e+04/2000) # FIXME

           if energy == 2458:
              scintEscale = (1/(1-lam))*(6.63328e+04*eff + ldarkrate)/2458
              scintE = NfinalPhotons/scintEscale # FIXME
              #ionizE = ThermalElectrons/(1.21347e+05/2458) # FIXME
              #noise = 50, threshold = 150.
              #ionizE = ThermalElectrons/(1.21892e+05/2458) # FIXME
              #noise = 50, threshold = 250.
              #ionizE = ThermalElectrons/(1.21803e+05/2458) # FIXME
              #noise = 100, threshold = 300.
              #ionizE = ThermalElectrons/(1.21814e+05/2458) # FIXME
              #noise = 100, threshold = 500.
              #ionizE = ThermalElectrons/(1.21678e+05/2458) # FIXME
              #noise = 200, threshold = 600.
              #ionizE = ThermalElectrons/(1.21602e+05/2458) # FIXME
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(1.21488e+05/2458) # FIXME
              #noise = 400, threshold = 1200.
              ionizE = ThermalElectrons/(1.21566e+05/2458) # FIXME

        if efield == 1000:
           if energy == 2458:
              scintEscale = (1/(1-lam))*(6.22660e+04*eff + ldarkrate)/2458 # SS
              scintE = NfinalPhotons/scintEscale
              #noise = 300, threshold = 900.
              #ionizE = ThermalElectrons/(1.25093e+05/2458) # FIXME 
              #noise = 200, threshold = 600.
              ionizE = ThermalElectrons/(1.25904e+05/2458) # FIXME 
#           if energy == 100:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(3.16161e+03*eff + ldarkrate)/100
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(3.82819e+03/100)
#           if energy == 200:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(5.83321e+03*eff + ldarkrate)/200
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(8.85683e+03/200)
#           if energy == 300:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(8.45643e+03*eff + ldarkrate)/300
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.39092e+04/300)
#           if energy == 400:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(1.10568e+04*eff + ldarkrate)/400
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.89976e+04/400)
#           if energy == 500:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(1.36315e+04*eff + ldarkrate)/500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(2.41030e+04/500) # FIXME
#           if energy == 1000:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(2.62989e+04*eff + ldarkrate)/1000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(4.98541e+04/1000)
#           if energy == 1500:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(3.87906e+04*eff + ldarkrate)/1500
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(7.57843e+04/1500)
#           if energy == 2000:
#              #noise = 200, threshold = 600.
#              scintEscale = (1/(1-lam))*(5.10642e+04*eff + ldarkrate)/2000
#              scintE = NfinalPhotons/scintEscale # FIXME
#              ionizE = ThermalElectrons/(1.01923e+05/2000)

        #################################################################
        ####   diffusion coefficient = 200 cm^2/s                     ###
        ################################################################# 
        #if efield == 65:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(1.05165e+05*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(8.30783e+04/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(8.29116e+04/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(8.27802e+04/2458)
        #      #noise = 300, threshold = 900. 
        #      ionizE = ThermalElectrons/(8.24435e+04/2458)

        #if efield == 130:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(9.36506e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(9.45792e+04/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(9.44320e+04/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(9.42496e+04/2458)
        #      #noise = 300, threshold = 900. 
        #      ionizE = ThermalElectrons/(9.39752e+04/2458)

        #if efield == 200:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(8.63334e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(1.01843e+05/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(1.01826e+05/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.01569e+05/2458)
        #      #noise = 300, threshold = 900. 
        #      ionizE = ThermalElectrons/(1.01363e+05/2458)

        #if efield == 300:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(7.96883e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(1.08523e+05/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(1.08398e+05/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.08198e+05/2458)
        #      #noise = 300, threshold = 900. 
        #      ionizE = ThermalElectrons/(1.08009e+05/2458)

        #if efield == 400:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(7.52264e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(1.13062e+05/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(1.12961e+05/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.12747e+05/2458)
        #      #noise = 300, threshold = 900. 
        #      ionizE = ThermalElectrons/(1.12464e+05/2458)

        #if efield == 500:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(7.19294e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(1.16317e+05/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(1.16201e+05/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.16007e+05/2458)
        #      #noise = 300, threshold = 900. 
        #      ionizE = ThermalElectrons/(1.15832e+05/2458)

        #if efield == 750:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(6.63328e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(1.21886e+05/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(1.21776e+05/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.21582e+05/2458)
        #      #noise = 300, threshold = 900. 
        #      ionizE = ThermalElectrons/(1.21396e+05/2458)

        #if efield == 1000:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(6.27453e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 50, threshold = 150. 
        #      #ionizE = ThermalElectrons/(1.25484e+05/2458)
        #      #noise = 100, threshold = 300. 
        #      #ionizE = ThermalElectrons/(1.25411e+05/2458)
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.24961e+05/2458)
        #      #noise = 300, threshold = 900. 
        #      #ionizE = ThermalElectrons/(1.24968e+05/2458)

        #################################################################
        ####   diffusion coefficient = 100 cm^2/s; 2458 MeV gamma; MS ###
        ################################################################# 
        #if efield == 65:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(1.11794e+05*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(7.41415e+04/2458)

        #if efield == 130:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(9.98370e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(8.63989e+04/2458)

        #if efield == 200:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(9.21227e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(9.43069e+04/2458)

        #if efield == 300:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(8.50221e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.01632e+05/2458)

        #if efield == 400:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(8.02017e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.06486e+05/2458)

        #if efield == 500:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(7.66148e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.10167e+05/2458)

        #if efield == 750:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(7.05530e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.16309e+05/2458)

        #if efield == 1000:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(6.66193e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      #ionizE = ThermalElectrons/(1.20296e+05/2458)

        #################################################################
        ####   diffusion coefficient = 100 cm^2/s; 2458 MeV gamma; SS ###
        ################################################################# 
        #if efield == 65:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(1.07711e+05*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(8.01976e+04/2458)

        #if efield == 130:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(9.60202e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(9.19019e+04/2458)

        #if efield == 200:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(8.83360e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(9.96204e+04/2458)

        #if efield == 300:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(8.17531e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(1.06341e+05/2458)

        #if efield == 400:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(7.71076e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(1.10901e+05/2458)

        #if efield == 500:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(7.37714e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(1.14357e+05/2458)

        #if efield == 750:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(6.79819e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(1.20142e+05/2458)

        #if efield == 1000:
        #   if energy == 2458:
        #      scintEscale = (1/(1-lam))*(6.41329e+04*eff + ldarkrate)/2458
        #      scintE = NfinalPhotons/scintEscale
        #      #noise = 200, threshold = 600. 
        #      ionizE = ThermalElectrons/(1.24016e+05/2458)

        scintE_tree[0] = NfinalPhotons
        ionizE_tree[0] = ThermalElectrons
        evtTree.Fill()

        hist['scintEraw'].Fill( scintE )
        hist['ionizEraw'].Fill( ionizE )

        # calibrate the energy
        #pI = [ -1/fittedTau, Zposition/velocity, 0. ]
        scintE = calibScintE(scintE, pS)
        ionizE = calibIonizE(ionizE, pI)
        hist['scintE'].Fill( scintE )
        hist['ionizE'].Fill( ionizE )
        hist['anti'].Fill( ionizE, scintE )

        EnergyList2D_ss.append( (ionizE, scintE) )

    ######################################################
#    hasElecNoise = '_NoElecNoise'
    hasElecNoise = ''

    global canvas
    global FileNames
    FileNames['SS'] = '%s/%s.pdf'%(dataDir, dataname)

    canvas = ROOT.TCanvas()
    for name in FileNames.itervalues(): canvas.Print(name + '[')

    # determine the optimal rotation angle and energy resolution
    ThetaToTry = [0.2 + 0.05*i for i in range(15)]
    #ThetaToTry = [0.0 + 0.01*i for i in range(100)]
    graph_ss = ROOT.TGraphErrors()
#    graph_ms = ROOT.TGraphErrors()

    BestTheta_ss = None
    BestRes_ss = None
#    BestTheta_ms = None
#    BestRes_ms = None

    global ssScintMean
    global ssIonizMean
    global msScintMean
    global msIonizMean
    ssScintResolution, ssScintMean = GetResolutionForTheta('SS Scint', EnergyList2D_ss, ROOT.TMath.Pi()/2)
    ssIonizResolution, ssIonizMean = GetResolutionForTheta('SS Ioniz', EnergyList2D_ss, 0.)
#    msScintResolution, msScintMean = GetResolutionForTheta('MS Scint', EnergyList2D_ms, ROOT.TMath.Pi()/2)
#    msIonizResolution, msIonizMean = GetResolutionForTheta('MS Ioniz', EnergyList2D_ms, 0.)

#    print ssIonizMean, ssScintMean
#    print msIonizMean, msScintMean
#    print ssIonizMean[0], ssScintMean[0]
#    print msIonizMean[0], msScintMean[0]

    for TestTheta in ThetaToTry:
        print "TestTheta ", TestTheta
        ssResolution, _ = GetResolutionForTheta('SS', EnergyList2D_ss, TestTheta)
#        msResolution, _ = GetResolutionForTheta('MS', EnergyList2D_ms, TestTheta)

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

    # ToDO:  Consider trying to pick a fit window more intelligently.
    ssPolynomialFunc = ROOT.TF1("ssFit", "[0] + [1]*(x-[2])**2", BestTheta_ss - 0.07, BestTheta_ss + 0.07)
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

    ssResolution, ssMean = GetResolutionForTheta('SS', EnergyList2D_ss, ssPolynomialFunc.GetParameter(2))
#    msResolution, msMean = GetResolutionForTheta('MS', EnergyList2D_ms, msPolynomialFunc.GetParameter(2))

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

   #######################################################
    c1 = ROOT.TCanvas('c1', 'c1', 1200, 500)
    c1.Divide(2,1)
    c1.cd(1)
    hist['scintEtrue'].Draw()
    hist['scintEtrue'].Fit('gaus', '+')
    c1.cd(2)
    #hist['ionizEtrue'].Draw()
    #hist['ionizEtrue'].Fit('gaus', '+')
    c1.Print('%s/%s_Etrue.eps'%(dataDir, dataname))

    c2 = ROOT.TCanvas('c2', 'c2', 1200, 500)
    c2.Divide(2,1)
    c2.cd(1)
    hist['scintEraw'].Draw()
    hist['scintEraw'].Fit('gaus', '+')
    c2.cd(2)
    hist['ionizEraw'].Draw()
    hist['ionizEraw'].Fit('gaus', '+')
    c2.Print('%s/%s_Eraw.eps'%(dataDir, dataname))

    c5 = ROOT.TCanvas('c5', 'c5', 1200, 500)
    c5.Divide(2,1)
    c5.cd(1)
    hist['NumScintGammas'].Draw()
    hist['NumScintGammas'].Fit('gaus', '+')
    c5.cd(2)
    hist['TECounter'].Draw()
    hist['TECounter'].Fit('gaus', '+')
    c5.Print('%s/%s_SignalRaw.eps'%(dataDir, dataname))

    c3 = ROOT.TCanvas('c3', 'c3', 600, 600)
    hist['anti'].Draw('COLZ')
    c3.Print('%s/%s_Anti.eps'%(dataDir,dataname))

    c4 = ROOT.TCanvas('c4', 'c4', 600, 500)
    hist['ionizE'].Draw()
    hist['scintE'].Draw('same')
    hist['scintE'].SetLineColor(2)
    c4.Print('%s/%s_calibE.eps'%(dataDir,dataname))

    output = ROOT.TFile('%s/%s.root'%(dataDir, dataname), 'recreate')
    for name in hist.iterkeys():
        hist[name].Write()
    evtTree.Write()
    output.Close()

    file = open('%s/%s.txt'%(dataDir, dataname), 'w')
    resNames = [ 'Theta_ss', 'PeakPos_ss', 'Resolution_ss',
                 'ssScintPeak', 'ssScintRes', 'ssIonizPeak', 'ssIonizRes' ]
    for name in resNames:
        value = result[name]
        if type(value) == type(1.0):
            file.write( str(float(value))+"\n")
        if type(value) == type(()):
            file.write( str(float(value[0]))+"  "+str(float(value[1]))+"\n"  )
    file.close()

    raw_input()
