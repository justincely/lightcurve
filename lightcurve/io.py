"""Library of I/O routines to get data into a LightCurve object.

"""

import matplotlib as mpl
#-- Don't render plots to screen
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import numpy as np

from .lightcurve import LightCurve
import io
from .cos import CosCurve
from .stis import StisCurve

__all__ = ['open']

#-------------------------------------------------------------------------------

def composite(filelist, output):
    """Creates a composite lightcurve from files in filelist and saves
    it to the save_loc.

    Parameters
    ----------
    filelist : list
        A list of full paths to the input files.
    output : string
        The path to the location in which the composite lightcurve is saved.
    """

    print("Creating composite lightcurve")

    wmin = 700
    wmax = 20000

    for filename in filelist:
        with pyfits.open(filename) as hdu:
            dq = hdu[1].data['DQ']
            wave = hdu[1].data['wavelength']
            xcorr = hdu[1].data['xcorr']
            ycorr = hdu[1].data['ycorr']
            sdqflags = hdu[1].header['SDQFLAGS']

            if (hdu[0].header['INSTRUME'] == "COS") and (hdu[0].header['DETECTOR'] == 'FUV'):
                other_file = [item for item in lightcurve.cos.get_both_filenames(filename) if not item == filename][0]
                if os.path.exists(other_file):
                    with pyfits.open(other_file) as hdu_2:
                        dq = np.hstack([dq, hdu_2[1].data['DQ']])
                        wave = np.hstack([wave, hdu_2[1].data['wavelength']])
                        xcorr = np.hstack([xcorr, hdu_2[1].data['xcorr']])
                        ycorr = np.hstack([ycorr, hdu_2[1].data['ycorr']])
                        sdqflags |= hdu_2[1].header['SDQFLAGS']

            index = np.where((np.logical_not(dq & sdqflags)) &
                             (wave > 500) &
                             (xcorr >=0) &
                             (ycorr >=0))

            if not len(index[0]):
                print('No Valid events')
                continue

            max_wave = wave[index].max()
            min_wave = wave[index].min()

            if max_wave < wmax:
                wmax = max_wave
            if min_wave > wmin:
                wmin = min_wave


    print('Using wavelength range of: {}-{}'.format(wmin, wmax))

    out_lc = LightCurve()

    for filename in filelist:
        new_lc = io.open(filename=filename,
                         step=1,
                         wlim=(wmin, wmax),
                         alt=None,
                         filter=True)

        out_lc = out_lc.concatenate(new_lc)

    out_lc.write(output, clobber=True)

#-------------------------------------------------------------------------------

def open(**kwargs):
    """ Open file into lightcurve

    filename must be supplied in kwargs

    Parameters
    ----------
    **kwargs : dict
        Additional arguements to be passed to lightcurve instantiations

    Returns
    -------
    LightCurve or subclass

    Raises
    ------
    IOError
        If filename is not supplied in **kwargs

    """

    if not 'filename' in kwargs:
        raise IOError('filename must be supplied')

    filetype = check_filetype(kwargs['filename'])

    if filetype == 'corrtag':
        return CosCurve(**kwargs)
    elif filetype == 'tag' or filetype == 'stis_corrtag':
        return StisCurve(**kwargs)
    elif filetype == 'lightcurve':
        return open_lightcurve(kwargs['filename'])
    else:
        raise IOError("Filetype not recognized: {}".format(filetype))

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
        filetype = 'corrtag'
    elif input_names == lightcurve_names:
        filetype = 'lightcurve'
    elif input_names == tag_names:
        filetype = 'tag'
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

    out_lc = LightCurve()

    hdu = pyfits.open(filename)

    out_lc.times = hdu[1].data['times']
    out_lc.gross = hdu[1].data['gross']
    out_lc.mjd = hdu[1].data['mjd']
    out_lc.flux = hdu[1].data['flux']
    out_lc.background = hdu[1].data['background']

    return out_obj

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
