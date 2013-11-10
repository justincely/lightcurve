"""
Tests for the LightCurve Object itsself.  

"""

from ..lightcurve import LightCurve
import numpy as np

#-------------------------------------------------------------------------------

def test_empty():
    """ Test that a new class instance contains empty arrays """

    obj = LightCurve()
    for item in [ obj.times, obj.mjd, obj.counts, obj.background]:
        assert np.array_equal( item, np.array( [] ) ), 'Arrays are not empty'
                
#-------------------------------------------------------------------------------

def test_net():
    """ Test the net calculation """

    obj = LightCurve()
    obj.counts = np.ones( 10 ) * 5
    obj.times = np.ones( 10 ) * 2

    assert np.array_equal( obj.net, np.ones( 10 ) * 2.5), \
        'Net not calculated right'
               
#-------------------------------------------------------------------------------

def test_error():
    """ Test the net calculation """

    obj = LightCurve()
    obj.counts = np.ones( 10 ) * 20
    obj.background = np.zeros( 10 )

    assert np.array_equal( obj.error, np.ones( 10 ) * np.sqrt( 20 )), \
        'Error not calculated right'

    obj.background = np.ones( 10 ) * 5

    assert np.array_equal( obj.error, np.ones( 10 ) * 5 ), \
        'Error not calculated right'

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    test_empty()
