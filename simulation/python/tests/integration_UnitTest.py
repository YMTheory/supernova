import unittest
import numpy as np

from task import visible_spectrum

class integration_UnitTest(unittest.TestCase):

    def test_dummy_integration(self):
        visible_spectrum.integrate()

