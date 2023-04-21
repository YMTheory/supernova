import unittest
import numpy as np

from xsec import NuE_XS
from xsec import NuP_XS
from xsec import IBD_XS

class xsecTest(unittest.TestCase):

    def test_creation(self):
        # --- Test class creation --- #
        try:
            self.nue_xsec = NuE_XS.NuE_XS()
        except:
            self.fail("NuE channel cross section class failed.")
        
        try:
            self.nup_xsec = NuP_XS.NuP_XS()
        except:
            self.fail("NuP channel cross section class failed.")
        
        try:
            self.ibd_xsec = IBD_XS.IBD_XS()
        except:
            self.fail("IBD channel cross section class failed.")

        return None

    # def test_dummy_diffxsec(self):
    #     # --- Test calculation of differential cross section --- #
    #     Enu = 10
    #     T   = np.arange(0, 10, 0.1)

    #     diffxs_NuP = np.zero_like(len(T))
    #     diffxs_NuE = np.zero_like(len(T))
    #     diffxs_IBD = np.zero_like(len(T))

    #     for tmpT in T:
    #         differential

