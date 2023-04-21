import numpy as np

from xsec import IBD_XS
from flux import Nakazato_2013_Model

from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.flavor_transformation import NoTransformation, AdiabaticMSW, ThreeFlavorDecoherence

def integrate():

    model_cls = Nakazato_2013_Model()
    model_cls.initialize_model()

    ibd = IBD_XS()
    
    def func_noDet(t, Ev):
        return ibd.totXS(Ev) * model_cls.get_arrvied_neutrino_flux()
    
    Emin, Emax = 0.8, 80
    E = np.linspace(Emin, Emax, 1000)
    tmin, tmax = -0.03, 0.07
    t = np.linspace(tmin, tmax, 100)

    result = np.trapz(func_noDetc, E, axis=0)

    print(result)
