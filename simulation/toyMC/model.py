import numpy as np
import uproot as up


def model(a, T, t0, t):
    if t0 < t < t0 + T:
        return 1 + a * np.sin(t)
    else:
        return 1



def generator(a, T, t0):

    t = np.arange(0, 1000, 0.1)
    Nevt = 10000
    Npoint = 1000
    y = np.zeros((Nevt, Npoint))

    for ievt in range(Nevt):
        for ipt in range(Npoint):
            t = np.random.uniform
            y[ievt, ipt] = model(a, T, t0, t)
    






