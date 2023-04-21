from flux.CCSNeModel import CCSNeModel

from astropy import units as u
from snewpy.neutrino import Flavor, MassHierarchy


class Nakazato_2013_Model(CCSNeModel):

    def __init__(self):
        super().__init__(model_name="Nakazato_2013")
        self.model_path = "/junofs/users/miaoyu/supernova/models/SNEWPY_models/Nakazato_2013/" 

        self.EoS    = "shen"
        self.z      = 0.004
        self.t_rev  = 100
        self.BH     = False
        self.s      = 20.0

    def initialize_model(self):
        
        from snewpy.models import Nakazato_2013

        path = "/junofs/users/miaoyu/supernova/models/SNEWPY_models/Nakazato_2013/"
        self.model_file = path + f"nakazato-{self.EoS}-z{self.z}-t_rev{self.t_rev:d}ms-s{self.s:.1f}.fits"
        self.model = Nakazato_2013(self.model_file)


    def get_flux_ET(self, t, Ev, flavor):
        """
        get luminosity (unit: [erg-1s-1]) at time t (unit: s) and neutrino energy E (unit: MeV) w/o flavor conversion
        """
        time = t * u.s
        E = Ev * u.MeV
        flu = self.model.get_initial_spectra(time, E)[flavor].to_value()
        
        return flu 


    def get_arrvied_neutrino_flux(self, t, Ev, distance, flavor, transfrom):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        time = t * u.s
        E = Ev * u.MeV
        return self.model.get_flux(time, E, distance, flavor_xform=transfrom)




