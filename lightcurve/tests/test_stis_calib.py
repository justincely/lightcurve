from ..stis import dqinit
from ..stis_calib import map_image

import numpy as np

def test_map_image():
    """Test the cython extension"""


    dq = np.random.random_integers(0, 10, 100)
    x, y = np.meshgrid(np.arange(10), np.arange(10))
    x = x.flatten()
    y = y.flatten()
    dq_im = dq.reshape((10, 10))


    #-- Test with 0s    
    assert np.array_equal(dq, map_image(dq_im, x, y)), "Mapping failed with ints"
