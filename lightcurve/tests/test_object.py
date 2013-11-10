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
    obj.times = np.ones( 10 ) * 2

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

    assert np.array_equal( obj.error, np.ones( 10 ) * np.sqrt( 20) ), \
        'Error not calculated right with a constant background'

#-------------------------------------------------------------------------------

def test_operations():
    """ test *, +, / """

    a = LightCurve()
    a.gross = np.ones( 2 ) * 10
    a.background = np.zeros( 2 )
    a.times = np.ones( 2 )
    a.mjd = a.times.copy()

    b = LightCurve()
    b.gross = np.ones( 3 ) * 20 
    b.background = np.zeros( 3 )
    b.times = np.ones( 3 )
    b.mjd = b.times.copy()

    assert np.array_equal( (a + 2).gross, np.ones( 2 ) * 10 + 2 ), \
        'Array addition by value not successful'

    assert np.array_equal( (a + b).gross, np.array( [10, 10, 20, 20, 20] ) ), \
        'Array concatenation not successful'

    assert np.array_equal( (a * 2).gross, np.ones( 2 ) * 10 * 2 ), \
        'Array multiplication not successful' 

    assert np.array_equal( (a / 2).gross, np.ones( 2 ) * 10 * .5 ), \
        'Array division not successful'

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    test_empty()
