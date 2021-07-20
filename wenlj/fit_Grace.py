#!/usr/bin/env python
# coding=utf-8

from ROOT import TGraphErrors, TF1, TCanvas, TMath, TFitResultPtr, TMatrixD, TGraph
import numpy as np
import matplotlib.pyplot as plt

def func(x, a0, aUV, aIR):
    return np.sqrt(a0 + aUV*x**2/(x**2-106.6**2) + aIR*x**2/(x**2-908.3**2) )

if __name__ == "__main__":


    wavelength = [361.2, 365, 406.3, 435.8, 475.3, 508.6, 546.1, 578, 643.9]
    #wavelength = np.array(wavelength) / 1000.  # um unit
    rindex = [1.2395, 1.2392, 1.2372, 1.2361, 1.2349, 1.2341, 1.2334, 1.2328, 1.2321]
    wavelength = np.array(wavelength)
    rindex = np.array(rindex)

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func, wavelength, rindex)

    print(popt)
    print(pcov)

    plotx = np.arange(300, 700, 1)
    plt.plot(wavelength, rindex,  "o", ms=5, label="Sinnock 1969")
    plt.plot(plotx, func(plotx, *popt), 'r-',
         label='fit: a0=%.3f, a1=%.4f, a2=%.6f' % tuple(popt))
    plt.xlabel('wavelength/nm')
    plt.ylabel('rindex')
    plt.show()


    """
    graph_data = TGraphErrors()
    graph_data2 = TGraphErrors() # n2 values
    #graph_data = TGraph()
    for idx, x, y in zip([i for i in range(9)], wavelength, rindex):
        graph_data.SetPoint(idx, x, y)
        graph_data.SetPointError(idx, 0, y*0.005)
        graph_data2.SetPoint(idx, x, y*y)
        graph_data2.SetPointError(idx, x, 2*y*0.005)


    # function form from Grace 2017
    func =  TF1("func", "([0]+[1]*x*x/(x*x-106.6*106.6)+[2]*x*x/(x*x-908.3*908.3))", 300, 700)
    #func =  TF1("func", "TMath::Sqrt([0]+[1]*x*x/(x*x-106.6*106.6)+[2]*x*x/(x*x-908.3*908.3))", 300, 700)
    func.SetParameter(0, 1.24)
    r = graph_data2.Fit(func, "RS")
    a0 = func.GetParameter(0)
    a1 = func.GetParameter(1)
    a2 = func.GetParameter(2)
    a0err = func.GetParError(0)
    a1err = func.GetParError(1)
    a2err = func.GetParError(2)

    cor = r.GetCorrelationMatrix();
    cov = r.GetCovarianceMatrix();

    #for i in range(3):
    #    for j in range(3):
    #        print(str(cor[i][j]), end=" ")
    #    print("\n")


    # function form from Sellmeier expression:
    #funcS = TF1("funcS", 
    #    "TMath::Sqrt(3/                                    \
    #        (1 - 1.2055/100*2/3*(35.49/(44.66/1000)) *     \
    #            ( [0]/(91.012-1/x/x)+                      \
    #              [1]/(87.892-1/x/x)+                      \
    #              [2]/(214.02-1/x/x) )                     \
    #        ) -2)", 0.3, 0.7)
    #funcS.FixParameter(0, 0.2075)
    #funcS.FixParameter(1, 0.0415)
    #funcS.FixParameter(2, 4.333)
    #graph_data.Fit(funcS, "RE")
    #a0 = funcS.GetParameter(0)
    #a1 = funcS.GetParameter(1)
    #a2 = funcS.GetParameter(2)
    #a0err = funcS.GetParError(0)
    #a1err = funcS.GetParError(1)
    #a2err = funcS.GetParError(2)

    drawx = np.arange(300, 700, 1)
    #drawx = np.arange(0.3, 0.7, 0.001)
    drawy = []
    for elem in drawx:
        drawy.append(np.sqrt(func.Eval(elem)) )

    plt.plot(drawx, drawy, "-", color='forestgreen', label="fitting w/ Grace's formula")
    plt.plot(wavelength, rindex, "o", color='darkviolet', label="Sinnock 1969")
    
    plt.text(0.55, 1.240, "a0: %.2f +- %.4f" %(a0, a0err))
    plt.text(0.55, 1.238, "a1: %.2f +- %.4f" %(a1, a1err))
    plt.text(0.55, 1.236, "a2: %.5f +- %.6f" %(a2, a2err))

    plt.xlabel('wavelength/nm')
    plt.ylabel('rindex')
    plt.legend()
    plt.show()
    """