import unittest

from snewpy.neutrino import Flavor, MassHierarchy
from snewpy.flavor_transformation import NoTransformation, AdiabaticMSW, ThreeFlavorDecoherence

from flux import Nakazato_2013_Model

class CCSNeModelTest(unittest.TestCase):

    def setUp(self):
        self.model = None

    def tearDown(self):
        self.model = None

    def _test_creation(self):
        # --- Test model creation --- #
        try:
            self.model = Nakazato_2013_Model.Nakazato_2013_Model()
            self.model.initialize_model()
        except:
            self.fail("Nakazato_2013_Model failed.")

    def test_get_initial_neutrino_flux_NU_E(self):
        self._test_creation()
        print(self.model.get_flux_ET_NU_E(0.01, 10))

    def test_get_arrived_neutrino_flux_noosc_NU_E(self):
        self._test_creation()
        flux = self.model.get_arrvied_neutrino_flux_noosc_NU_E(0.01, 10)
        print(f"Time = 0.01s, neutrino energy = 10 MeV, flux from 10 kpc CCSNe without flavour conversion = {flux}.")
        
    def test_get_arrived_neutrino_flux_NO_NU_E(self):
        self._test_creation()
        flux = self.model.get_arrvied_neutrino_flux_NO_NU_E(0.01, 10)
        print( f"Time = 0.01s, neutrino energy = 10 MeV, flux from 10 kpc CCSNe with NO flavour conversion = {flux}. ")

    def test_get_arrived_neutrino_flux_IO_NU_E(self):
        self._test_creation()
        flux = self.model.get_arrvied_neutrino_flux_IO_NU_E(0.01, 10)
        print( f"Time = 0.01s, neutrino energy = 10 MeV, flux from 10 kpc CCSNe with IO flavour conversion = {flux} [cm^-2 erg^-1 s^-1]. ")


