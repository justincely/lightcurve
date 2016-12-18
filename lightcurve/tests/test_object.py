"""
Tests for the LightCurve Object itsself.

"""
from astropy.table import Table
import numpy as np

from .. import io

#-------------------------------------------------------------------------------

class test_derived_columns:
    def setUp(self):
        self.data = {'times': np.arange(10),
                     'mjd': np.arange(10),
                     'bins': np.ones(10),
                     'gross': np.ones(10)*5,
                     'background': np.ones(10)*1,
                     'flux': np.ones(10)*1}

        self.lc = io.read(self.data)

    def tearDown(self):
        pass

    def test_counts(self):
        assert np.array_equal(self.lc['counts'], np.ones(10)*4), \
            'Counts not calculated right'

    def test_net(self):
        assert np.array_equal(self.lc['net'], np.ones(10)*4), \
            'Net not calculated right'

    def test_error(self):
        assert np.array_equal(self.lc['error'], np.ones(10) * np.sqrt(6)), \
            'Error not calculated right with no background'

#-------------------------------------------------------------------------------
