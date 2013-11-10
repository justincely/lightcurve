from ..lightcurve import LightCurve
import numpy as np

def test_empty():
    """ Test that a new class instance contains empty arrays """
    obj = LightCurve()
    for item in [ obj.times, obj.mjd, obj.counts, obj.background]:
        assert np.array_equal( item, np.array( [] ) ), 'Arrays are not empty'
                               
if __name__ == "__main__":
    test_empty()
