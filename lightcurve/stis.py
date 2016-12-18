""" Utility functions for extracting STIS spectral data into lightcurves

"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import os
import numpy as np
import scipy
from scipy.interpolate import interp1d
from datetime import datetime
from numba import jit
import astropy
from astropy.io import fits as fits

from .utils import expand_refname, enlarge
from .cos import extract_index, calc_npixels
from .version import version as  __version__

__all__ = ['extract',
           'stis_corrtag',
           'map_image',
           'epsilon',
           'dqinit']

#-------------------------------------------------------------------------------

def extract(filename, **kwargs):
    """ Extract lightcurve from STIS dataset

    This is the main driver of the lightcuve extracion, and definitely
    needs some better documentation.

    Parameters
    ----------
    filename : str
        name of FITS file to extract from
    **kwargs : dict
        arbitrary keyword arguements for tailored extraction

    Kwarg parameters
    ----------------
    verbosity : int, default=0
        Verbosity level for print output
    step : int, default=1
        timestep in seconds for output Lightcurve
    wlim : tuple

    Returns
    -------
    data, meta : Astropy table, dict
        Table with extracted data and dictionary of metadata pairs
    """

    verbosity = kwargs.get('verbosity', 0)
    step = kwargs.get('step', 1)
    wlim = kwargs.get('wlim', (2, 10000))
    xlim = kwargs.get('xlim', (0, 2048))
    ylim = kwargs.get('ylim', (0, 2048))
    filter_airglow = kwargs.get('filter_airglow', True)

    hdu = fits.open(filename)
    if fits.getval(filename, 'OBSTYPE') == 'IMAGING':
        print("Imaging observation found, resetting limits.")
        xlim = (0, 2048)
        ylim = (0, 2048)
        wlim = (-1, 1)

    input_headers = {'a':{}}

    for i, ext in enumerate(hdu):
        try:
            input_headers['a'][i] = ext.header._cards
        except AttributeError:
            pass

    meta = {'source': filename,
            'instrument' : 'STIS',
            'headers': input_headers,
            'stepsize': step,
            'wlim': wlim,
            'xlim': xlim,
            'ylim': ylim}

    if verbosity:
        print('#-----------------------------------------#')
        print('Running LightCurve extraction for STIS Data')
        print('#-----------------------------------------#')
        print()
        print('Extracting from: {}'.format(filename))
        print('With arguments:')
        print('step : {}'.format(step))
        print('wlim : {}'.format(wlim))
        print('xlim : {}'.format(xlim))
        print('ylim : {}'.format(ylim))
        print('filter_airglow : {}'.format(filter_airglow))
        print()

    SECOND_PER_MJD = 1.15741e-5

    time = hdu[1].data['time']
    #time = np.array([round(val, 3) for val in hdu[1].data['time']]).astype(np.float64)

    if not len(time):
        end = 0
        start = 0
    else:
        end = min(time.max(), hdu[1].header['EXPTIME'])
        start = time.min()

    all_steps = np.arange(start, end+step, step)

    if all_steps[-1] > end:
        truncate = True
    else:
        truncate = False

    ### need to fix sdqflags
    index = extract_index(hdu,
                          xlim[0], xlim[1],
                          ylim[0], ylim[1],
                          wlim[0], wlim[1],
                          hdu[1].header['SDQFLAGS'],
                          filter_airglow=filter_airglow)

    if verbosity:
        print("#{} events".format(len(index)))

    gross = np.histogram(hdu['events'].data['time'][index],
                         all_steps,
                         weights=hdu['events'].data['epsilon'][index])[0]

    n_pixels = calc_npixels(hdu, index, xlim)

    response_array = get_fluxes(hdu, index).mean()
    weights = hdu['events'].data['epsilon'][index] / step / response_array
    print("WARNING: The flux is a lie")
    flux =  np.histogram(hdu['events'].data['time'][index], all_steps, weights=weights)[0] / n_pixels

    print("WARNING: The background is a lie")
    background = np.zeros(gross.shape)
    mjd = hdu[1].header['EXPSTART'] + np.array(all_steps[:-1]) * SECOND_PER_MJD
    bins = np.ones(len(gross)) * step
    times = all_steps[:-1]

    if truncate:
        if verbosity:
            print('Truncating the last event bin')

        gross = gross[:-1]
        flux = flux[:-1]
        background = background[:-1]
        mjd = mjd[:-1]
        bins = bins[:-1]
        times = times[:-1]

    data = {'dataset': np.ones(times.shape),
            'times': times,
            'mjd': mjd,
            'bins': bins,
            'gross': gross,
            'background': background,
            'flux': flux}

    if verbosity:
        print('Finished extraction for {}'.format(filename))
        print()

    return data, meta

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
    if header0['OBSTYPE'] == 'IMAGING':
        wave_data = np.zeros(n_events)

    elif not os.path.exists(x1d_filename):
        print("Could not find associated extracted spectrum {}".format(x1d_filename))
        wave_data = np.ones(n_events) * hdu[0].header['CENTRWV']

    elif header0['OPT_ELEM'].startswith('E'):
        #-- Put in average wavelength for now
        with fits.open(x1d_filename) as x1d:
            #--initialize everything to 0 wavelength
            wave_data = np.zeros(n_events)

            offset = header1['shifta2']
            for order in x1d[1].data:
                spec_center = order['a2center'] + offset
                spec_height = order['extrsize']

                index = np.where((ycorr_data < spec_center+spec_height) &
                                 (ycorr_data > spec_center-spec_height))[0]

                int_pix = integerize_pixels(xcorr_data[index])
                wave_data[index] = order['wavelength'][int_pix]

    else:
        #-- Grab wavelengths from the x1d file
        int_pix = integerize_pixels(xcorr_data)
        with fits.open(x1d_filename) as x1d:
            wave_data = x1d[1].data['wavelength'][0][int_pix]

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

def integerize_pixels(xcoords):

    int_pix = np.round(xcoords).astype(np.int32)
    int_pix = np.where(int_pix < 0, 0, int_pix)
    int_pix = np.where(int_pix > 2047, 2047, int_pix)
    int_pix //= 2

    return int_pix

#-------------------------------------------------------------------------------

@jit
def map_image(image, xcoords, ycoords, default=0):
    n_coord = len(xcoords)
    out_vals = np.zeros(n_coord)

    for i in range(n_coord):
        x = xcoords[i]
        y = ycoords[i]

        if x < 0 or x >= 2048:
            val = default
        elif y < 0 or y >= 2048:
            val = default
        else:
            val = image[y, x]

        out_vals[i] = val

    return out_vals

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

            if not os.path.exists(reffile):
                print("{} not found, correction not performed".format(reffile))
                return np.ones(len(hdu[1].data))

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

        if not os.path.exists(reffile):
            print("{} not found, correction not performed".format(reffile))
            return np.zeros(len(hdu[1].data))

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
    grid = scipy.sparse.coo_matrix((weights, xyi), shape=(nx, ny)).toarray().T

    return grid, np.linspace(xmin, xmax, nx), np.linspace(ymin, ymax, ny)

#-------------------------------------------------------------------------------

def get_fluxes(hdu, index):
    """ Return response curve and wavelength range for specified mode
    """

    fluxtab = hdu[0].header.get('FLUXTAB', '')

    fluxfile = expand_refname(fluxtab)

    if not len(index):
        return np.ones(hdu['events'].data['time'].shape)[index]

    if (not fluxfile) or (not os.path.exists(fluxfile)):
        print(' WARNING: Fluxfile not available %s,' % fluxfile )
        print(' using unity flux calibration instead.')
        return np.ones(hdu['events'].data['time'].shape)[index]


    flux_hdu = fits.open(fluxfile)
    setting_index = np.where((flux_hdu[1].data['SPORDER'] == hdu[0].header['segment']) &
                             (flux_hdu[1].data['OPT_ELEM'] == hdu[0].header['opt_elem']) &
                             (flux_hdu[1].data['CENWAVE'] == hdu[0].header['cenwave']))[0]

    if len(setting_index) == 0:
        print('No row in fluxtab found for this dataset, no FLUXCAL performed')
        return np.ones(hdu['events'].data['time'].shape)[index]
    elif len(setting_index) > 1:
        raise ValueError('Too many rows found: {}'.format(len(setting_index)))


    resp_wave = flux_hdu[1].data[setting_index]['WAVELENGTH'].flatten()
    response = flux_hdu[1].data[setting_index]['SENSITIVITY'].flatten()

    data_max = hdu['events'].data['wavelength'][index].max()
    data_min = hdu['events'].data['wavelength'][index].min()

    if data_max > resp_wave.max():
        print("Expanding minumum response curve by {}".format(data_max - resp_wave.max()))
        resp_wave[np.argmax(resp_wave)] = data_max

    if data_min < resp_wave.min():
        print("Expanding minimum response curve by {}".format(data_min - resp_wave.min()))
        resp_wave[np.argmin(resp_wave)] = data_min

    interp_func = interp1d(resp_wave, response)
    all_resp = interp_func(hdu['events'].data['wavelength'][index])

    return all_resp

#-------------------------------------------------------------------------------
