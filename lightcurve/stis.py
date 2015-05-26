"""
Holder of Space Telescope Imaging Spectrograph (STIS) classes and utilities

"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import os
import numpy as np
import scipy
from scipy.interpolate import interp1d
from datetime import datetime
import astropy
from astropy.io import fits as fits

from .utils import expand_refname, enlarge
from .stis_cal import map_image
from .cos import extract_index
from .version import version as  __version__

#-------------------------------------------------------------------------------

def extract(filename, **kwargs):
    """ Loop over HDUs and extract the lightcurve

    This is the main driver of the lightcuve extracion, and definitely
    needs some better documentation.

    """

    step = kwargs.get('step', 1)
    wlim = kwargs.get('wlim', (2, 10000))
    xlim = kwargs.get('xlim', (0, 2048))
    ylim = kwargs.get('ylim', (0, 2048))
    filter_airglow = kwargs.get('filter_airglow', True)

    hdu = fits.open(filename)

    meta = {'source': filename,
            'hdus': {'A': hdu},
            'stepsize': step,
            'wlim': wlim,
            'xlim': xlim,
            'ylim': ylim}

    print("Using wlim: {}".format(wlim))
    print("Using xlim: {}".format(xlim))
    print("Using ylim: {}".format(ylim))
    print("      step: {}".format(step))

    SECOND_PER_MJD = 1.15741e-5


    time = hdu[1].data['time']
    #time = np.array([round(val, 3) for val in hdu[1].data['time']]).astype(np.float64)

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

    ### need to fix sdqflags
    index = extract_index(hdu,
                          xlim[0], xlim[1],
                          ylim[0], ylim[1],
                          wlim[0], wlim[1],
                          0,
                          filter_airglow=filter_airglow)

    print("#{} events".format(len(index)))

    gross = np.histogram(hdu['events'].data['time'][index],
                          all_steps,
                          weights=hdu['events'].data['epsilon'][index])[0]

    flux = np.zeros(gross.shape)
    background = np.zeros(gross.shape)
    mjd = hdu[1].header['EXPSTART'] + np.array(all_steps[:-1]) * SECOND_PER_MJD
    bins = np.ones(len(gross)) * step
    times = all_steps[:-1]

    if truncate:
        gross = gross[:-1]
        flux = flux[:-1]
        background = background[:-1]
        mjd = mjd[:-1]
        bins = bins[:-1]
        times = times[:-1]

    data = [times, mjd, bins, gross, background, flux]
    columns = ('times', 'mjd' ,'bins', 'gross', 'background', 'flux')

    return data, columns, meta

#-------------------------------------------------------------------------------

def stis_corrtag(tagfile, clean=True):
    """Create a COS-like corrtag file for STIS data

    Parameters
    ----------
    tagfile, str
        input STIS time-tag data file

    """
    print("Creating STIS corrtag for {}".format(tagfile))
    if '_tag.fits' in tagfile:
        x1d_filename = tagfile.replace('_tag.fits', '_x1d.fits')
    elif '_corrtag.fits' in tagfile:
        x1d_filename = tagfile.replace('_corrtag.fits', '_x1d.fits')

    with fits.open(tagfile, 'readonly') as hdu:
        n_events = len(hdu[1].data['TIME'])

        time_data = hdu[1].data['TIME']

        if 'DETAXIS1' in hdu[1].data.names:
            #-- Y data is same in raw and corrected coordinates
            rawx_data = hdu[1].data['DETAXIS1']
            rawy_data = hdu[1].data['AXIS2']
            xcorr_data = hdu[1].data['AXIS1']
            ycorr_data = hdu[1].data['AXIS2']
        else:
            rawx_data = hdu[1].data['RAWX']
            rawy_data = hdu[1].data['RAWY']
            xcorr_data = hdu[1].data['XCORR'].astype(np.int32)
            ycorr_data = hdu[1].data['YCORR'].astype(np.int32)

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
        with fits.open(x1d_filename) as x1d:
            #--initialize everything to 0 wavelength
            wave_data = np.zeros(n_events)

            #for order in x1d[1].data:
            #    spec_center = order['a2center'] - order['shifta2']
            #    spec_height = order['extrsize']

    else:
        #-- Grab wavelengths from the x1d file
        int_pix = np.array(map(int, map(round, xcorr_data))).astype(np.int32)
        int_pix = np.where(int_pix < 0, 0, int_pix)
        int_pix = np.where(int_pix > 2047, 2047, int_pix)
        int_pix //= 2

        with fits.open(x1d_filename) as x1d:
            wave_data = x1d[1].data['wavelength'][0][int_pix]
            print(wave_data.mean())

    #-- Writeout corrtag file
    hdu_out = fits.HDUList(fits.PrimaryHDU())

    hdu_out[0].header['GEN_DATE'] = (str(datetime.now()), 'Creation Date')
    hdu_out[0].header['LC_VER'] = (__version__, 'lightcurve version used')
    hdu_out[0].header['AP_VER'] = (astropy.__version__, 'Astropy version used')
    hdu_out[0].header['NP_VER'] = (np.__version__, 'Numpy version used')
    hdu_out[0].header['SP_VER'] = (scipy.__version__, 'Scipy version used')

    hdu_out[0].header.extend(header0, end=True)

    time_col = fits.Column('time', 'D', 'second', array=time_data)
    rawx_col = fits.Column('rawx', 'D', 'pixel', array=rawx_data)
    rawy_col = fits.Column('rawy', 'D', 'pixel', array=rawy_data)
    xcorr_col = fits.Column('xcorr', 'D', 'pixel', array=xcorr_data)
    ycorr_col = fits.Column('ycorr', 'D', 'pixel', array=ycorr_data)
    epsilon_col = fits.Column('epsilon', 'D', 'factor', array=eps_data)
    wave_col = fits.Column('wavelength', 'D', 'angstrom', array=wave_data)
    dq_col = fits.Column('dq', 'I', 'DQ', array=dq_data)

    tab = fits.new_table([time_col,
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
    if clean and '_tag.fits' in tagfile:
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

    with fits.open(tagfile) as hdu:

        epsilon_out = np.ones(hdu[1].data['time'].shape)

        #-- Flatfield correction
        for ref_flat in ['PFLTFILE', 'LFLTFILE']:
            reffile = expand_refname(hdu[0].header[ref_flat])
            print('FLATFIELD CORRECTION {}: {}'.format(ref_flat, reffile))

            with fits.open(reffile) as image_hdu:
                image = image_hdu[1].data

                if not image.shape == (2048, 2048):
                    x_factor = 2048 // image.shape[1]
                    y_factor = 2048 // image.shape[0]

                    print('Enlarging by {},{}'.format(x_factor, y_factor))
                    image = enlarge(image, x_factor, y_factor)

                #--indexing is 1 off
                if 'AXIS1' in hdu[1].data.names:
                    epsilon_out *= map_image(image,
                                             hdu[1].data['AXIS1'] - 1,
                                             hdu[1].data['AXIS2'] - 1)
                else:
                    epsilon_out *= map_image(image,
                                             hdu[1].data['XCORR'].astype(np.int32) - 1,
                                             hdu[1].data['YCORR'].astype(np.int32) - 1)

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

    with fits.open(tagfile) as hdu:

        reffile = expand_refname(hdu[0].header['BPIXTAB'])
        print('BPIXTAB used {}'.format(reffile))

        with fits.open(reffile) as bpix:
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
        if 'AXIS1' in hdu[1].data.names:
            dq_out = map_image(dq_im,
                               hdu[1].data['AXIS1'] - 1,
                               hdu[1].data['AXIS2'] - 1)
        else:
            dq_out = map_image(dq_im,
                               hdu[1].data['XCORR'].astype(np.int32) - 1,
                               hdu[1].data['YCORR'].astype(np.int32) - 1)

    return dq_out

#-------------------------------------------------------------------------------

def crazy_histogram2d(x, y, bins=(2048, 2048)):
    """  Bin events to a 2d image.

    This is faster than the numpy version.

    http://stackoverflow.com/questions/8805601/efficiently-create-2d-histograms-from-large-datasets

    Parameters
    ----------
    x : array-like
        x-coordinates to bin
    y : array-like
        y-coordinates to bin
    bins : tuple, optional
        bins of output image

    Returns
    -------
    grid, x_bins, y_bins : (np.ndarray, np.ndarray, np.ndarray)
        binned image, x bins, and y bins
    """

    try:
        nx, ny = bins
    except TypeError:
        nx = ny = bins

    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    dx = (xmax - xmin) / (nx - 1.0)
    dy = (ymax - ymin) / (ny - 1.0)

    weights = np.ones(x.size)

    # Basically, this is just doing what np.digitize does with one less copy
    xyi = np.vstack((x,y)).T
    xyi = np.floor(xyi, xyi).T

    # Now, we'll exploit a sparse coo_matrix to build the 2D histogram...
    grid = scipy.sparse.coo_matrix((weights, xyi), shape=(nx, ny)).toarray()

    return grid, np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)
