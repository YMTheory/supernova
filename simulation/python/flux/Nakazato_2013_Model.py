from flux.CCSNeModel import CCSNeModel

from astropy import units as u
from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.flavor_transformation import NoTransformation, AdiabaticMSW, ThreeFlavorDecoherence

import warnings
warnings.filterwarnings("ignore")


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


    def get_flux_ET_NU_E(self, t, Ev):
        """
        get luminosity (unit: [erg-1s-1]) at time t (unit: s) and neutrino energy E (unit: MeV) w/o flavor conversion
        """
        time = t * u.s
        E = Ev * u.MeV
        flu = self.model.get_initial_spectra(time, E)[Flavor.NU_E].to_value()
        
        return flu 

    def get_flux_ET_NU_E_BAR(self, t, Ev):
        """
        get luminosity (unit: [erg-1s-1]) at time t (unit: s) and neutrino energy E (unit: MeV) w/o flavor conversion
        """
        time = t * u.s
        E = Ev * u.MeV
        flu = self.model.get_initial_spectra(time, E)[Flavor.NU_E_BAR].to_value()
        
        return flu 

    def get_flux_ET_NU_X(self, t, Ev):
        """
        get luminosity (unit: [erg-1s-1]) at time t (unit: s) and neutrino energy E (unit: MeV) w/o flavor conversion
        """
        time = t * u.s
        E = Ev * u.MeV
        flu = self.model.get_initial_spectra(time, E)[Flavor.NU_X].to_value()
        
        return flu 

    def get_flux_ET_NU_X_BAR(self, t, Ev):
        """
        get luminosity (unit: [erg-1s-1]) at time t (unit: s) and neutrino energy E (unit: MeV) w/o flavor conversion
        """
        time = t * u.s
        E = Ev * u.MeV
        flu = self.model.get_initial_spectra(time, E)[Flavor.NU_X_BAR].to_value()
        
        return flu 


    def get_arrvied_neutrino_flux_noosc_NU_E(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        transfrom = NoTransformation()
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=transfrom)[Flavor.NU_E]

    def get_arrvied_neutrino_flux_noosc_NU_E_BAR(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        transfrom = NoTransformation()
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=transfrom)[Flavor.NU_E_BAR]

    def get_arrvied_neutrino_flux_noosc_NU_X(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        transfrom = NoTransformation()
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=transfrom)[Flavor.NU_X]

    def get_arrvied_neutrino_flux_noosc_NU_X_BAR(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        transfrom = NoTransformation()
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=transfrom)[Flavor.NU_X_BAR]



    def get_arrvied_neutrino_flux_NO_NU_E(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.NORMAL)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_E]

    def get_arrvied_neutrino_flux_NO_NU_E_BAR(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.NORMAL)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_E_BAR]

    def get_arrvied_neutrino_flux_NO_NU_X(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.NORMAL)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_X]

    def get_arrvied_neutrino_flux_NO_NU_X_BAR(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.NORMAL)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_X_BAR]


    def get_arrvied_neutrino_flux_IO_NU_E(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.INVERTED)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_E]

    def get_arrvied_neutrino_flux_IO_NU_E_BAR(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.INVERTED)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_E_BAR]

    def get_arrvied_neutrino_flux_IO_NU_X(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.INVERTED)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_X]

    def get_arrvied_neutrino_flux_IO_NU_X_BAR(self, t, Ev):
        """
        get the neutrino number (unit: /(cm^2*erg*s)) at time t and neutrino energy Ev from distance away.
        """
        nmo_trans = AdiabaticMSW(mh=MassHierarchy.INVERTED)
        time = t * u.s
        E = Ev * u.MeV
        distance = 10.  # unit: kpc
        return self.model.get_flux(time, E, distance, flavor_xform=nmo_trans)[Flavor.NU_X_BAR]

