import numpy as np
import uproot as up

def get_statROI_from1Dpdf(filename:str, tmin=-0.02, tmax=0.02) -> float:
    """
    get statistics in [tmin, tmax] time window from 1D pdf.
    """
    f = up.open(filename)
    h = f["h1"]
    x = h.axis().centers()
    y = h.values()
    xwidth = x[1] - x[0]
    nevt = 0
    for i, j in zip(x, y):
        if i >= tmin and i <= tmax:
            nevt += j * xwidth
    return nevt

def get_statROI_from2Dpdf(filename:str, tmin=-0.02, tmax=0.02, Emin=0.10, Emax=5):
    """
    get statistics in ROI in [tmin, tmax] and [Emin, Emax] range from 2D pdf.
    :param filename:
    :param tmin:
    :param tmax:
    :param Emin:
    :param Emax:
    :return:
    """
    pass

def get_statROI_from1Ddata(filename: str):
    """
    get statistics mean value from 1D data.
    :param filename: input data file
    :return: mean value of event in data file
    """
    f = up.open(filename)
    arr = f["binned"]["TbinConts"].array()
    Nevt = [len(subarr) for subarr in arr]
    return np.mean(Nevt)


def get_statROI_from2Ddata(filename:str):
    """
    get statistics mean value from 2D dataset.
    :param filename: input filename.
    :return: mean value of event in data file.
    """
    f = up.open(filename)
    arr = f["binned"]["TbinConts"].array()
    Nevt = [len(subarr) for subarr in arr]
    return np.mean(Nevt)