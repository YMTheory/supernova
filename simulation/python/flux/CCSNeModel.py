from astropy import units as u

from abc import ABC, abstractmethod

class CCSNeModel(ABC):

    def __init__(self, model_name):
        self.model_name = model_name
        self.model_path = "./"
        self.model = None

    @abstractmethod
    def initialize_model(self):
        pass

    @abstractmethod
    def get_flux_ET(self, t, Ev, flavor):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux(self, t, Ev, distance, flavor, transfrom):
        return 0.










