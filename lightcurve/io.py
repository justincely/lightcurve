"""Library of I/O routines to get data into a LightCurve object.

"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import os
import scipy
import numpy as np
from datetime import datetime
import astropy
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.time import Time

from .version import version as __version__
from .cos import extract as extract_cos
from .cos import get_both_filenames
from .stis import extract as extract_stis
from .utils import is_uniq

__all__ = ['check_filetype',
           'read',
           'composite',
           'prepare_header']

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

    corrtag_names = set(['TIME',
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
                         'PHA'])

    lightcurve_names = set(['TIMES',
                            'MJD',
                            'GROSS',
                            'COUNTS',
                            'NET',
                            'FLUX',
                            'FLUX_ERROR',
                            'BACKGROUND',
                            'ERROR'])

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

    with fits.open(filename) as hdu:
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

def read(source=None, **kwargs):
    verbosity = kwargs.get('verbosity', 0)

    if isinstance(source, dict):
        data = source
        meta = {}

    elif source == None:
        if verbosity: print("Initializing empty table.")

        data = {'dataset': np.array([]),
                'times': [],
                'mjd': [],
                'bins': [],
                'gross': [],
                'background': [],
                'flux': [],
                'counts': [],
                'error': [],
                'flux_error': [],
                'net': [],
                'signal_to_noise': []}

        meta = {'filename': None}

    else:
        filetype = check_filetype(source)
        if verbosity: print("Found {} as input filetype:".format(filetype))

        if filetype == 'lightcurve':
            return Table.read(source)
        elif filetype == 'cos_corrtag':
            data, meta = extract_cos(source, **kwargs)
        elif filetype == 'stis_tag' or filetype == 'stis_corrtag':
            data, meta = extract_stis(source, **kwargs)
        else:
            raise IOError("Filetype not recognized: {}".format(filetype))

        meta['outname'] = source[:9] + '_curve.fits'
        meta['filetype'] = filetype

        meta['GEN_DATE'] = (str(datetime.now()), 'Creation Date')
        meta['LC_VER'] = (__version__, 'lightcurve version used')
        meta['AP_VER'] = (astropy.__version__, 'Astropy version used')
        meta['NP_VER'] = (np.__version__, 'Numpy version used')
        meta['SP_VER'] = (scipy.__version__, 'Scipy version used')

        meta['expstart'] = data['mjd'].min()
        meta['expend'] = data['mjd'].max()
        meta['exptime'] = data['bins'].sum()
        meta['stepsize'] = (meta['stepsize'], 'Bin size (seconds)')

        meta['WMIN'] = (meta['wlim'][0], 'Minimum wavelength extracted')
        meta['WMAX'] = (meta['wlim'][1], 'Maximum wavelength extracted')
        meta['XMIN'] = (meta['xlim'][0], 'Minimum x-coordinate extracted')
        meta['XMAX'] = (meta['xlim'][1], 'Maximum x-coordinate extracted')

        if meta['instrument'] == 'COS':
            if 'ylima' in meta:
                meta['YMIN_A'] = (meta['ylima'][0], 'Minimum y-coordinate extracted from FUVA')
                meta['YMAX_A'] = (meta['ylima'][1], 'Maximum y-coordinate extracted from FUVA')
            if 'ylimb' in meta:
                meta['YMIN_B'] = (meta['ylimb'][0], 'Minimum y-coordinate extracted from FUVB')
                meta['YMAX_B'] = (meta['ylimb'][1], 'Maximum y-coordinate extracted from FUVB')
        elif meta['instrument'] == 'STIS':
            meta['YMIN'] = (meta['ylim'][0], 'Minimum y-coordinate extracted')
            meta['YMAX'] = (meta['ylim'][1], 'Maximum y-coordinate extracted')
        else:
            raise ValueError("{} not a recognized instrument".format(meta['instrument']))

    data['counts'] = data['gross'] - data['background']
    data['error'] = np.sqrt(data['gross'] + data['background'])
    data['flux_error'] = np.sqrt(data['gross'] + data['background']) / (data['gross'] - data['background']) * data['flux']
    data['net'] = (data['gross'] - data['background']) / data['bins']
    data['signal_to_noise'] = data['gross'] / np.sqrt(data['gross'] + data['background'])

    return Table(data, meta=meta)

#-------------------------------------------------------------------------------

def quicklook(filename):
    """ Quick plotting function for extracted lightcurves
    """

    with fits.open(filename) as hdu:

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(hdu[1].data['times'], hdu[1].data['gross'], 'o')

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Gross Counts')

        fig.suptitle(filename)
        fig.savefig(filename.replace('.fits', '.png'))
        plt.close(fig)

#-------------------------------------------------------------------------------

def composite(filelist, output, trim=True, **kwargs):
    """Creates a composite lightcurve from files in filelist and saves
    it to the save_loc.

    Parameters
    ----------
    filelist : list
        A list of full paths to the input files.
    output : string
        The path to the location in which the composite lightcurve is saved.
    trim : bool, opt
        Trim wavelengths to common ranges for all files
    """

    print("Creating composite lightcurve from:")
    print("\n".join(filelist))

    wmin = 912
    wmax = 9000

    for filename in filelist:
        with fits.open(filename) as hdu:
            dq = hdu[1].data['DQ']
            wave = hdu[1].data['wavelength']
            xcorr = hdu[1].data['xcorr']
            ycorr = hdu[1].data['ycorr']
            sdqflags = hdu[1].header['SDQFLAGS']

            if (hdu[0].header['INSTRUME'] == "COS") and (hdu[0].header['DETECTOR'] == 'FUV'):
                #-- if using COS FUV, these are the good limits
                wmin = 912
                wmax = 1800
                other_file = [item for item in get_both_filenames(filename) if not item == filename][0]
                if os.path.exists(other_file):
                    with fits.open(other_file) as hdu_2:
                        dq = np.hstack([dq, hdu_2[1].data['DQ']])
                        wave = np.hstack([wave, hdu_2[1].data['wavelength']])
                        xcorr = np.hstack([xcorr, hdu_2[1].data['xcorr']])
                        ycorr = np.hstack([ycorr, hdu_2[1].data['ycorr']])
                        sdqflags |= hdu_2[1].header['SDQFLAGS']

            index = np.where((np.logical_not(dq & sdqflags)) &
                             (wave > 500) &
                             (xcorr >= 0) &
                             (ycorr >= 0))

            if not len(index[0]):
                print('No Valid events')
                continue

            wmax = min(wmax, wave[index].max())
            wmin = max(wmin, wave[index].min())
            print(wmin, '-->', wmax)


    if trim:
        print('Using wavelength range of: {}-{}'.format(wmin, wmax))
        kwargs['wlim'] = (wmin, wmax)

    all_lc = []
    for filename in filelist:
        print(filename)
        tmp_lc = read(filename, **kwargs)
        if np.any(tmp_lc['gross'] == 0):
            continue
        else:
            all_lc.append(tmp_lc)

    for i, lc in enumerate(all_lc):
        if i == 0:
            out_lc = lc
        else:
            lc['dataset'] = i+1
            out_lc = vstack([out_lc, lc], metadata_conflicts='warn')

    if output.endswith('.fits') or output.endswith('.fits.gz'):
        for key in out_lc.meta.keys():
            if len(key) > 8:
                print("Deleting key {} from output header:".format(key))
                del out_lc.meta[key]

    out_lc.write(output)

    prepare_header(output, filelist)

#-------------------------------------------------------------------------------

def prepare_header(filename, filelist):
    """Prepare headers with MAST requirements"""
    telescop = set()
    instrume = set()
    detector = set()
    filter = set()

    for i, exposure in enumerate(filelist):
        with fits.open(exposure) as hdu:
            telescop.add(hdu[0].header['TELESCOP'])
            instrume.add(hdu[0].header['INSTRUME'])
            detector.add(hdu[0].header['DETECTOR'])
            filter.add(hdu[0].header['OPT_ELEM'])

            if i == 0:
                ra_targ = hdu[0].header['ra_targ']
                dec_targ = hdu[0].header['dec_targ']
                equinox = hdu[0].header['equinox']
                targname = hdu[0].header['targname']
                tardescr = hdu[0].header.get('TARDESCR', '')
                tardesc2 = hdu[0].header.get('TARDESC2', '')

    with fits.open(filename, mode='update') as hdu:
        #-- HSLP keywords
        hdu[0].header['PROPOSID'] = 13902
        hdu[0].header['HLSPLEAD'] = 'Justin C. Ely'
        hdu[0].header['PR_INV_L'] = 'Ely'
        hdu[0].header['PR_INV_F'] = 'Justin'
        hdu[0].header['PR_INV_M'] = 'Charles'
        hdu[0].header['HLSPNAME'] = 'The Lightcurve Legacy of COS and STIS'
        hdu[0].header['HLSPACRN'] = 'LLOCS'
        hdu[0].header['CITATION'] = ''

        hdu[0].header['RA_TARG'] = ra_targ
        hdu[0].header['DEC_TARG'] = dec_targ
        hdu[0].header['EQUINOX'] = equinox
        hdu[0].header['TARGNAME'] = targname
        hdu[0].header['TARDESCR'] = tardescr
        hdu[0].header['TARDESC2'] = tardesc2
        hdu[0].header['EXTNAME'] = 'PRIMARY'
        hdu[1].header['EXTNAME'] = 'LIGHTCURVE'

        uniq, value = is_uniq(telescop)
        hdu[0].header['telescop'] = value

        uniq, value = is_uniq(instrume)
        hdu[0].header['instrume'] = value
        if not uniq:
            for i, val in enumerate(list(instrume)):
                hdu[0].header['instru{:0>2}'.format(1)] = val

        uniq, value = is_uniq(detector)
        hdu[0].header['detector'] = value
        if not uniq:
            for i, val in enumerate(list(detector)):
                hdu[0].header['detect{:0>2}'.format(1)] = val


        uniq, value = is_uniq(filter)
        hdu[0].header['filter'] = value
        if not uniq:
            for i, val in enumerate(list(filter)):
                hdu[0].header['filter{:0>2}'.format(1)] = val

        hdu[0].header['DATE-OBS'] = Time(hdu[1].data['MJD'].min(), format='mjd').iso
        hdu[0].header['EXPSTART'] = hdu[1].data['mjd'].min()
        hdu[0].header['EXPEND'] = hdu[1].data['mjd'].max()
        hdu[0].header['EXPTIME'] = hdu[1].data['bins'].sum()

#-------------------------------------------------------------------------------
