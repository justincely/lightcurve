"""Library of I/O routines to get data into a LightCurve object.

"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import matplotlib as mpl
#-- Don't render plots to screen
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import numpy as np
import os

#-------------------------------------------------------------------------------

def check_filetype(filename):
    """Determine the type of data being input.

    File type is determined by the culumns in the first data extension.

    Parameters
    ----------
    filename : str
        name of the input file

    Returns
    -------
    filetype : str
        determined type of file

    """

    corrtag_names = set( ['TIME',
                          'RAWX',
                          'RAWY',
                          'XCORR',
                          'YCORR',
                          'XDOPP',
                          'XFULL',
                          'YFULL',
                          'WAVELENGTH',
                          'EPSILON',
                          'DQ',
                          'PHA'] )

    lightcurve_names = set( ['TIMES',
                             'MJD',
                             'GROSS',
                             'COUNTS',
                             'NET',
                             'FLUX',
                             'FLUX_ERROR',
                             'BACKGROUND',
                             'ERROR'] )

    tag_names = set(['TIME',
                     'AXIS1',
                     'AXIS2',
                     'DETAXIS1'])

    stis_corrtag_names = set(['TIME',
                              'RAWX',
                              'RAWY',
                              'XCORR',
                              'YCORR',
                              'WAVELENGTH',
                              'EPSILON',
                              'DQ'])

    with pyfits.open(filename) as hdu:
        input_names = set([item.upper() for
                           item in hdu[1].data.names])

    if input_names == corrtag_names:
        filetype = 'cos_corrtag'
    elif input_names == lightcurve_names:
        filetype = 'lightcurve'
    elif input_names == tag_names:
        filetype = 'stis_tag'
    elif input_names == stis_corrtag_names:
        filetype = 'stis_corrtag'
    else:
        filetype = None

    return filetype

#-------------------------------------------------------------------------------

def open_lightcurve(filename):
    """ Read lightcurve from fits file back into base object

    Parameters
    ----------
    filename : str
        Filename of FITS lightcurve to open

    Returns
    -------
    out_lc : LightCurve
        LightCurve instantiation containing data from file

    """

    with pyfits.open(filename) as hdu:
        columns = hdu[1].data.names
        data = [hdu[1].data[name] for name in columns]
        meta = {}

    return data, columns, meta

#-------------------------------------------------------------------------------

def quicklook(filename):
    """ Quick plotting function for extracted lightcurves
    """

    with pyfits.open(filename) as hdu:

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(hdu[1].data['times'], hdu[1].data['gross'], 'o')

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Gross Counts')

        fig.suptitle(filename)
        fig.savefig(filename.replace('.fits', '.png'))
        plt.close(fig)

#-------------------------------------------------------------------------------
