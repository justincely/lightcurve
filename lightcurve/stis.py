"""
Holder of Space Telescope Imaging Spectrograph (STIS) classes and utilities

"""

from __future__ import print_function

import os
import numpy as np
import scipy
from scipy.interpolate import interp1d
from datetime import datetime
import astropy
from astropy.io import fits as pyfits

from .lightcurve import LightCurve
from .utils import expand_refname, enlarge
from .stis_cal import calculate_epsilon, map_dq_image
from .version import version as  __version__

#-------------------------------------------------------------------------------

class StisCurve(LightCurve):
    """
    Subclass for STIS specific routines and abilities

    """

    def __init__(self, **kwargs):
        """ Initialize and extract lightcurve from input corrtag
        """

        if 'step' in kwargs: step = kwargs['step']
        else: step = 1

        if 'filename' in kwargs:
            self.extract(kwargs['filename'], step=step)


    def __str__(self):
        """Prettier representation of object instanct"""

        return "Lightcurve Object of %s" % ( ','.join( self.input_list ) )


    def extract(self, filename, step=1):
        """ Loop over HDUs and extract the lightcurve

        This is the main driver of the lightcuve extracion, and definitely
        needs some better documentation.

        """

        SECOND_PER_MJD = 1.15741e-5

        hdu = pyfits.open(filename)
        self.hdu = hdu

        time = hdu[1].data['time']
        #time = np.array([round(val, 3) for val in self.hdu[1].data['time']]).astype(np.float64)

        if not len(time):
            end = 0
        else:
            end = min(time.max(),
                      hdu[1].header['EXPTIME'] )


        all_steps = np.arange(0, end+step, step)

        if all_steps[-1] > end:
            truncate = True
        else:
            truncate = False

        gross = np.histogram(time, all_steps)[0]

        self.gross = gross
        self.flux = np.zeros(gross.shape)
        self.background = np.zeros(gross.shape)
        self.mjd = hdu[1].header['EXPSTART'] + np.array(all_steps[:-1]) * SECOND_PER_MJD
        self.bins = np.ones(len(gross)) * step
        self.times = all_steps[:-1]

        if truncate:
            self.gross = self.gross[:-1]
            self.flux = self.flux[:-1]
            self.background = self.background[:-1]
            self.mjd = self.mjd[:-1]
            self.bins = self.bins[:-1]
            self.times = self.times[:-1]

#-------------------------------------------------------------------------------

def stis_corrtag(tagfile, clean=True):
    """Create a COS-like corrtag file for STIS data

    Parameters
    ----------
    tagfile, str
        input STIS time-tag data file

    """
    print("Creating STIS corrtag for {}".format(tagfile))
    x1d_filename = tagfile.replace('_tag.fits', '_x1d.fits')

    with pyfits.open(tagfile, 'readonly') as hdu:
        n_events = len(hdu[1].data['TIME'])

        time_data = hdu[1].data['TIME']

        #-- Y data is same in raw and corrected coordinates
        rawx_data = hdu[1].data['DETAXIS1']
        rawy_data = hdu[1].data['AXIS2']
        xcorr_data = hdu[1].data['AXIS1']
        ycorr_data = hdu[1].data['AXIS2']

        #-- copy over the primary
        header0 = hdu[0].header.copy()
        header1 = hdu[1].header.copy()

    eps_data = epsilon(tagfile)
    dq_data = dqinit(tagfile)
    if not os.path.exists(x1d_filename):
        print("Could not find associated extracted spectrum {}".format(x1d_filename))
        wave_data = np.ones(n_events) * hdu[0].header['CENTRWV']
    elif header0['OPT_ELEM'].startswith('E'):
        #-- Put in average wavelength for now
        with pyfits.open(x1d_filename) as x1d:
            wave_data = np.ones(n_events) * x1d[1].data['wavelength'].mean()
    else:
        #-- Grab wavelengths from the x1d file
        int_pix = np.array(map(int, map(round, xcorr_data))).astype(np.int32)
        int_pix = np.where(int_pix < 0, 0, int_pix)
        int_pix = np.where(int_pix > 1023, 1023, int_pix)
        print(int_pix)
        with pyfits.open(x1d_filename) as x1d:
            wave_data = x1d[1].data['wavelength'][0][int_pix]
            print(wave_data.mean())

    #-- Writeout corrtag file
    hdu_out = pyfits.HDUList(pyfits.PrimaryHDU())

    hdu_out[0].header['GEN_DATE'] = (str(datetime.now()), 'Creation Date')
    hdu_out[0].header['LC_VER'] = (__version__, 'lightcurve version used')
    hdu_out[0].header['AP_VER'] = (astropy.__version__, 'Astropy version used')
    hdu_out[0].header['NP_VER'] = (np.__version__, 'Numpy version used')
    hdu_out[0].header['SP_VER'] = (scipy.__version__, 'Scipy version used')

    hdu_out[0].header.extend(header0, end=True)

    time_col = pyfits.Column('time', 'D', 'second', array=time_data)
    rawx_col = pyfits.Column('rawx', 'D', 'pixel', array=rawx_data)
    rawy_col = pyfits.Column('rawy', 'D', 'pixel', array=rawy_data)
    xcorr_col = pyfits.Column('xcorr', 'D', 'pixel', array=xcorr_data)
    ycorr_col = pyfits.Column('ycorr', 'D', 'pixel', array=ycorr_data)
    epsilon_col = pyfits.Column('epsilon', 'D', 'factor', array=eps_data)
    wave_col = pyfits.Column('wavelength', 'D', 'angstrom', array=wave_data)
    dq_col = pyfits.Column('dq', 'I', 'DQ', array=dq_data)

    tab = pyfits.new_table([time_col,
                            rawx_col,
                            rawy_col,
                            xcorr_col,
                            ycorr_col,
                            epsilon_col,
                            wave_col,
                            dq_col])
    hdu_out.append(tab)

    hdu_out[1].header.extend(header1, end=True)

    hdu_out.writeto(tagfile.replace('_tag.fits', '_corrtag.fits'), clobber=True)
    if clean:
        print("Removing input tagfile")
        os.remove(tagfile)
    
#-------------------------------------------------------------------------------

def epsilon(tagfile):
    """Compute the total epsilon factor for each event

    Compute the flatfield correction from the P-flat and L-flat reference files
    (PFLTFILE and LFLTFILE respectively).

    Parameters
    ----------
    tagfile, str
        input STIS time-tag data file

    Returns
    -------
    epsilon, np.ndarray
        array of epsilons

    """

    print("Calculating Epsilon")

    with pyfits.open(tagfile) as hdu:

        epsilon_out = np.ones(hdu[1].data['time'].shape)

        #-- Flatfield correction
        for ref_flat in ['PFLTFILE', 'LFLTFILE']:
            reffile = expand_refname(hdu[0].header[ref_flat])
            print('FLATFIELD CORRECTION {}: {}'.format(ref_flat, reffile))

            with pyfits.open(reffile) as image_hdu:
                image = image_hdu[1].data

                if not image.shape == (2048, 2048):
                    x_factor = 2048 / image.shape[1]
                    y_factor = 2048 / image.shape[0]

                    print('Enlarging by {},{}'.format(x_factor, y_factor))
                    image = enlarge(image, x_factor, y_factor)

                #--indexing is 1 off
                epsilon_out *= calculate_epsilon(image,
                                                 hdu[1].data['AXIS1'] - 1,
                                                 hdu[1].data['AXIS2'] - 1)

    return epsilon_out

#-------------------------------------------------------------------------------

def dqinit(tagfile):
    """Compute the data quality information for each pixel from the BPIXTAB.

    Parameters
    ----------
    tagfile, str
        input STIS time-tag data file

    Returns
    -------
    dq, np.ndarray
        array of bitwise dq flags

    """

    print("Calculating DQ")

    with pyfits.open(tagfile) as hdu:

        reffile = expand_refname(hdu[0].header['BPIXTAB'])
        print('BPIXTAB used {}'.format(reffile))

        with pyfits.open(reffile) as bpix:
            #-- Mama bpix regions are in lo-res pixels
            dq_im = np.zeros((1024, 1024)).astype(np.int32)

            #-- indexing is 1 off
            for line in bpix[1].data:
                lx = line['PIX1'] - 1
                ly = line['PIX2'] - 1
                dpix = line['LENGTH']
                axis = line['AXIS']
                flag = line['VALUE']

                if axis == 1:
                    #-- Along X
                    dq_im[ly, lx:lx + dpix] |= flag
                elif axis == 2:
                    #-- Along Y
                    dq_im[ly:ly+dpix, lx] |= flag

            #-- Needs to be expanded into Hi-res
            dq_im = enlarge(dq_im, 2, 2)


        #-- Map to the events
        #-- indexing is 1 off
        dq_out = map_dq_image(dq_im,
                              hdu[1].data['AXIS1'] - 1,
                              hdu[1].data['AXIS2'] - 1)

    return dq_out
