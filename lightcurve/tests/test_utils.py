""" Test utility functions in util module

"""

import os
import numpy as np

from ..utils import expand_refname, is_uniq, enlarge

from ..cos import expand_refname, extract_index

#-------------------------------------------------------------------------------

def test_expand_refname():
    """  Test that i can accurately find the header reference file """

    if 'lref' in os.environ: del os.environ['lref']
    assert expand_refname('lref$something_tds.fits') == 'something_tds.fits', \
        "Without lref, didn't return None"

    os.environ['lref'] = '/grp/hst/cdbs/lref/'
    full_refname = os.path.join( os.environ['lref'], 'something_tds.fits')
    assert expand_refname( 'lref$something_tds.fits' ) == full_refname, \
        "With lref, didn't find the file"

    assert expand_refname( '' ) == '', "didn't return blank on blank"

#-------------------------------------------------------------------------------

def test_is_uniq():
    empty_set = set()
    uniq_set = set([1])
    non_uniq_set = set([1, 2, 3])

    assert is_uniq(empty_set) == (True, ''), "Empty set failed"

    assert is_uniq(uniq_set) == (True, 1), "Unique set failed"

    assert is_uniq(non_uniq_set) == (False, 'MULTI'), "heterogeneous set failed"

#-------------------------------------------------------------------------------

def test_enlarge():
    small_array = np.ones((4, 4))
    big_array = np.ones((8, 8))

    assert np.array_equal(enlarge(small_array), big_array), "array enlarging failed"

#-------------------------------------------------------------------------------
