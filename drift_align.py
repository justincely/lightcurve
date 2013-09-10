#! /user/bin/env python
"""
Check COS data for uncorrected drift during a single exposure
"""

import pyfits
import numpy as np
from scipy.ndimage.filters import convolve

from ttag_funcs import ttag_image
from astroraf.spectra.spectools import fft_correlate

import sys ### change to argparse at some point.

def split_corrtag( filename, n_spectra=2 ):
    """ Split and extract corrtag data
    """

    hdu = pyfits.open( filename )

    start = hdu['events'].data['time'].min()
    end = hdu['events'].data['time'].max()
    step = (end - start) / n_spectra

    spectra_list = []
    
    while start < end-1:
        print 'Binning from %3.2f to %3.2f'% (start, start + step)

        image = ttag_image( filename, xtype='XFULL', ytype='YFULL',
                            times=(start, start+step) )
        
        spectra_list.append( np.sum(image, axis=0) )
        
        start += step

    return spectra_list


def check_for_drift( filename, n_spectra):
    import pylab
    spectra_list = split_corrtag( filename, n_spectra )

    for i in range(len( spectra_list )):
        spectra_list[i] = convolve( spectra_list[i], np.ones(12)/12, mode='mirror' )
    
    template_spectra = spectra_list[0]

    pylab.figure()
    pylab.plot( template_spectra )

    found_shifts = []
    for spectrum in spectra_list[1:]:
        shift = fft_correlate( template_spectra, spectrum )
        found_shifts.append( shift )
        
        pylab.plot( spectrum )
        

    print found_shifts
    if np.any( found_shifts ):
        raw_input()
    pylab.close()


if __name__ == "__main__":
    main( sys.argv[1], sys.argv[2] )
