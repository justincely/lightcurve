"""
Tests for the LightCurve Object itsself.  

"""

from ..lightcurve import LightCurve
import numpy as np

#-------------------------------------------------------------------------------

def test_empty():
    """ Test that a new class instance contains empty arrays """

    obj = LightCurve()
    for item in [ obj.times, obj.mjd, obj.gross, obj.counts, obj.background]:
        assert np.array_equal( item, np.array( [] ) ), 'Arrays are not empty'
                
#-------------------------------------------------------------------------------

def test_counts():
    """ Test the count calculation """

    obj = LightCurve()
    obj.gross = np.ones( 10 ) * 5
    obj.background = np.ones( 10 ) * 1 

    assert np.array_equal( obj.counts, np.ones( 10 ) * 4), \
        'Counts not calculated right'

#-------------------------------------------------------------------------------

def test_net():
    """ Test the net calculation """

    obj = LightCurve()
    obj.gross = np.ones( 10 ) * 5
    obj.background = np.zeros( 10 )
    obj.bins = np.ones( 10 ) * 2

    assert np.array_equal( obj.net, np.ones( 10 ) * 2.5), \
        'Net not calculated right'
               
#-------------------------------------------------------------------------------

def test_error():
    """ Test the error calculation """

    obj = LightCurve()
    obj.gross = np.ones( 10 ) * 20
    obj.background = np.zeros( 10 )

    assert np.array_equal( obj.error, np.ones( 10 ) * np.sqrt( 20 )), \
        'Error not calculated right with no background'

    obj.background = np.ones( 10 ) * 5
    assert np.array_equal( obj.error, np.ones( 10 ) * np.sqrt( 20 + 5 ) ), \
        'Error not calculated right with a constant background'

#-------------------------------------------------------------------------------

def test_operations():
    """ test *, +, / """

    a = LightCurve()
    a.gross = np.ones(2) * 10
    a.flux = np.ones(2)
    a.background = np.zeros(2)
    a.bins = np.ones(2)
    a.times = np.arange(2)
    a.mjd = a.bins.copy()

    b = LightCurve()
    b.gross = np.ones(3) * 20
    b.flux = np.ones(3) 
    b.background = np.zeros(3)
    b.bins = np.ones(3)
    b.times = np.arange(3)
    b.mjd = b.bins.copy()

    assert np.array_equal( (a + b).gross, np.array( [10, 10, 20, 20, 20] ) ), \
        'Array concatenation not successful'

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    test_empty()
