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
    def get_flux_ET_NU_E(self, t, Ev):
        return 0.

    @abstractmethod
    def get_flux_ET_NU_E_BAR(self, t, Ev):
        return 0.

    @abstractmethod
    def get_flux_ET_NU_X(self, t, Ev):
        return 0.

    @abstractmethod
    def get_flux_ET_NU_X_BAR(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_noosc_NU_E(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_noosc_NU_E_BAR(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_noosc_NU_X(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_noosc_NU_X_BAR(self, t, Ev):
        return 0.


    @abstractmethod
    def get_arrvied_neutrino_flux_NO_NU_E(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_NO_NU_E_BAR(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_NO_NU_X(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_NO_NU_X_BAR(self, t, Ev):
        return 0.


    @abstractmethod
    def get_arrvied_neutrino_flux_IO_NU_E(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_IO_NU_E_BAR(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_IO_NU_X(self, t, Ev):
        return 0.

    @abstractmethod
    def get_arrvied_neutrino_flux_IO_NU_X_BAR(self, t, Ev):
        return 0.











