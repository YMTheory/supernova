import numpy as np

from xsec import IBD_XS
from flux import Nakazato_2013_Model
from detector import toyDetector

from snewpy.neutrino import Flavor

def integrate():

    model_cls = Nakazato_2013_Model.Nakazato_2013_Model()
    model_cls.initialize_model()

    ibd = IBD_XS.IBD_XS()

    det = toyDetector.toyDetector(20, 0.12, 0.88)
    
    def func_noDet_NO(t, Ev):
        return ibd.totXS(Ev, Flavor.NU_E_BAR) * model_cls.get_arrvied_neutrino_flux_NO_NU_E_BAR(t, Ev) * det.nH
    ufunc_noDet_NO = np.vectorize(func_noDet_NO)

    tmp_val = ufunc_noDet_NO(0.02, 10)

    print(f"Middle check variable = {tmp_val:.2f}.")
    
    #Emin, Emax = 0.8, 80
    #E = np.linspace(Emin, Emax, 1000)
    #tmin, tmax = -0.03, 0.07
    #t = np.linspace(tmin, tmax, 100)

    #f = ufunc_noDet_NO(t, E)
    #print(f)
    #result = np.trapz(func_noDet_NO, E, axis=0)

