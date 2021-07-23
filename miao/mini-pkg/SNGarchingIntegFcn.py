import numpy as np
import ROOT
import sys
import re

datapath = "/junofs/users/miaoyu/supernova/wenlj/simulation/data/Garching/"

def readFile(filename):
    t, lum, aveE, aveE2 = [], [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = re.split('\ +',line)
            t.append(float(data[1]))
            lum.append(float(data[2]))
            aveE.append(float(data[3]))
            aveE2.append(float(data[4]))
    return t, lum, aveE, aveE2

grLuminosity = [ROOT.TGraph() for i in range(3)]
grAlpha = [ROOT.TGraph() for i in range(3)]
grAverageE = [ROOT.TGraph() for i in range(3)]
timeMin = [0 for i in range(3)]
timeMax = [0 for i in range(3)]

def readFluxGraph(imode):
    if imode < 70000 or imode > 99999:
        print("Wrong mode name %d, does not belong to Garching models!!" %imode)
        sys.exit(-1)

    prefix = str(imode) + "/timedata/neutrino_signal_nu"
    t_nue, lumin_nue, averagE_nue, averagE2_nue = readFile(datapath+prefix+"_e")
    t_antinue, lumin_antinue, averagE_antinue, averagE2_antinue = readFile(datapath+prefix+"bar_e")
    t_nux, lumin_nux, averagE_nux, averagE2_nux = readFile(datapath+prefix+"_x")

    # time limits:
    timeMin[0] = t_nue[0]
    timeMin[1] = t_antinue[0]
    timeMin[2] = t_nux[0]
    timeMax[0] = t_nue[-1]
    timeMax[1] = t_antinue[-1]
    timeMax[2] = t_nux[-1]

    for i in range(len(t_nue)):
        grLuminosity[0].SetPoint(i, t_nue[i], lumin_nue[i])
        grAlpha[0].SetPoint(i, t_nue[i], 1./(averagE2_nue[i]/averagE_nue[i]**2-1)-1)
        grAverageE[0].SetPoint(i, t_nue[i], averagE_nue[i])

    for i in range(len(t_antinue)):
        grLuminosity[1].SetPoint(i, t_antinue[i], lumin_antinue[i])
        grAlpha[1].SetPoint(i, t_antinue[i], 1./(averagE2_antinue[i]/averagE_antinue[i]**2-1)-1)
        grAverageE[1].SetPoint(i, t_antinue[i], averagE_antinue[i])

    for i in range(len(t_nux)):
        grLuminosity[2].SetPoint(i, t_nux[i], lumin_nux[i])
        grAlpha[2].SetPoint(i, t_nux[i], 1./(averagE2_nux[i]/averagE_nux[i]**2-1)-1)
        grAverageE[2].SetPoint(i, t_nux[i], averagE_nux[i])

    print("Load Garching flux mode %d successfully !!!"%imode)


from scipy.special import gamma
def GarchingFluxFcn(time, Ev, tp):
    luminosity, A, alpha = 0, 0, 0
    if tp < 2:
        luminosity = grLuminosity[tp].Eval(time)
        A = grAverageE[tp].Eval(time)
        alpha = grAlpha[tp].Eval(time)
    if tp >= 2:
        luminosity = grLuminosity[2].Eval(time)
        A = grAverageE[2].Eval(time)
        alpha = grAlpha[2].Eval(time)

    index = 6.24151e56
    flux = (luminosity / A) * (np.power(Ev, alpha) /  gamma(1+alpha)) * np.power((alpha+1)/A, alpha+1) * np.exp(-(alpha+1)*Ev/A)

    return index * flux


def getEventAtTime(time, Ev, tp):
    luminosity, A, alpha = 0, 0, 0
    if tp < 2:
        if time < timeMin[tp]:
            return 0
        if time > timeMax[tp]:
            return 0
        luminosity = grLuminosity[tp].Eval(time)
        A = grAverageE[tp].Eval(time)
        alpha = grAlpha[tp].Eval(time)
    if tp >= 2:
        if time < timeMin[2]:
            return 0
        if time > timeMax[2]:
            return 0
        luminosity = grLuminosity[2].Eval(time)
        A = grAverageE[2].Eval(time)
        alpha = grAlpha[2].Eval(time)

    index = 6.24151e56
    flux = (luminosity / A) * (np.power(Ev, alpha) /  gamma(1+alpha)) * np.power((alpha+1)/A, alpha+1) * np.exp(-(alpha+1)*Ev/A)

    return index * flux
    

def getAverageET(time, tp):
    aveE = 0
    if tp < 2:
        if timeMin[tp] < time < timeMax[tp]:
            aveE = grAverageE[tp].Eval(time)
        else:
            aveE = 0
    else:
        if timeMin[2] < time < timeMax[2]:
            aveE = grAverageE[2].Eval(time)
        else:
            aveE = 0

    return aveE


def getNumT(time, tp):
    num = 0
    if tp<2:
        if time < timeMin[tp]:
            return 0
        if time < timeMax[tp]:
            return 0
        num = grLuminosity[tp].Eval(time)/(grAverageE[tp].Eval(time))
    if tp >=2:
        if time < timeMin[2]:
            return 0
        if time < timeMax[2]:
            return 0
        num = grLuminosity[2].Eval(time)/(grAverageE[2].Eval(time))

    index = 6.24151e56
    return index * num



