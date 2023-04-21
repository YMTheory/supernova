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

    def test_get_initial_neutrino_flux(self):
        self._test_creation()
        for flavor in Flavor:
            self.model.get_flux_ET(10, 0.01, flavor)
            #print(self.model.get_flux_ET(10, 0.01, flavor))

    def test_get_arrived_neutrino_flux(self):
        self._test_creation()

        xform_nmo = AdiabaticMSW()
        for flavor in Flavor:
            print(self.model.get_arrvied_neutrino_flux(10, 0.01, 10, flavor, xform_nmo)[flavor])
        
