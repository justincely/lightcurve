"""
Tests for the LightCurve Object itsself.

"""
from astropy.table import Table
import numpy as np

from ..lightcurve import LightCurve

#-------------------------------------------------------------------------------

def test_empty():
    """ Test that a new class instance contains empty arrays """

    obj = LightCurve()
    for item in [obj['times'], obj['mjd'], obj['gross'], obj['background'], obj['flux']]:
        assert len(item) == 0, 'Arrays are not empty'

#-------------------------------------------------------------------------------

def test_counts():
    """ Test the count calculation """

    columns = ('times',
              'mjd',
              'bins',
              'gross',
              'background',
              'flux')

    data = [np.arange(10),
            np.arange(10),
            np.ones(10),
            np.ones(10)*5,
            np.ones(10)*1,
            np.ones(10)*1]

    test_tab = Table(data=data, names=columns, meta={})

    lc = LightCurve(test_tab)

    assert np.array_equal(lc.counts, np.ones(10)*4), \
        'Counts not calculated right'

#-------------------------------------------------------------------------------

def test_net():
    """ Test the net calculation """

    columns = ('times',
               'mjd',
               'bins',
               'gross',
               'background',
               'flux')

    data = [np.arange(10),
            np.arange(10),
            np.ones(10)*2,
            np.ones(10)*5,
            np.ones(10)*0,
            np.ones(10)*1]

    test_tab = Table(data=data, names=columns, meta={})

    lc = LightCurve(test_tab)

    assert np.array_equal(lc.net, np.ones(10)*2.5), \
        'Net not calculated right'

#-------------------------------------------------------------------------------

def test_error():
    """ Test the error calculation """

    columns = ('times',
               'mjd',
               'bins',
               'gross',
               'background',
               'flux')

    data = [np.arange(10),
            np.arange(10),
            np.ones(10)*2,
            np.ones(10)*20,
            np.ones(10)*0,
            np.ones(10)*1]

    test_tab = Table(data=data, names=columns, meta={})

    lc = LightCurve(test_tab)

    assert np.array_equal(lc.error, np.ones(10) * np.sqrt(20)), \
        'Error not calculated right with no background'

    lc['background'] = 5
    assert np.array_equal(lc.error, np.ones(10) * np.sqrt(20 + 5)), \
        'Error not calculated right with a constant background'

#-------------------------------------------------------------------------------

def test_operations():
    """ test *, +, / """


    columns = ('times',
               'mjd',
               'bins',
               'gross',
               'background',
               'flux')

    data = [np.arange(2),
            np.arange(2),
            np.ones(2),
            np.ones(2)*10,
            np.zeros(2),
            np.ones(2)]
    test_tab = Table(data=data, names=columns, meta={})
    a = LightCurve(test_tab)


    data = [np.arange(3),
            np.arange(3),
            np.ones(3),
            np.ones(3)*20,
            np.zeros(3),
            np.ones(3)]
    test_tab = Table(data=data, names=columns, meta={})
    b = LightCurve(test_tab)

    c = a.concatenate(b)
    assert np.array_equal(c['gross'], np.array([10, 10, 20, 20, 20])), \
        'Array concatenation not successful'

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    test_empty()
