"""
Contains the basic LightCurve object which is inherited by the subclasses.

"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

__all__ = ['LightCurve']

import astropy
import scipy
import numpy as np
from datetime import datetime
from astropy.io import fits
from astropy.table import Table

from .version import version as __version__
from .io import check_filetype, open_lightcurve
from .cos import extract as extract_cos
from .stis import extract as extract_stis

#-------------------------------------------------------------------------------

class LightCurve(Table):
    """
    Returns
    -------
    lightcurve object

    """

    def __init__(self, filename=None, **kwargs):
        """ Instantiate an empty LightCurve object.

        LightCurves must contain:
            gross
            bins (timestep)
            times
            mjd
            background
            flux

        All values are currently set to empty arrays.

        """

        if filename is None:
            print("Initializing emtpy LightCurve")
            data = [[], [], [], [], []]
            columns = ('times',
                            'mjd',
                            'bins',
                            'gross',
                            'background')
            meta = {'filename': None}

        else:
            filetype = check_filetype(filename)

            if filetype == 'cos_corrtag':
                data, columns, meta = extract_cos(filename, **kwargs)
            elif filetype == 'stis_tag' or filetype == 'stis_corrtag':
                data, columns, meta = extract_stis(filename, **kwargs)
            elif filetype == 'lightcurve':
                data, columns, meta = open_lightcurve(filename)
            else:
                raise IOError("Filetype not recognized: {}".format(filetype))


        super(LightCurve, self).__init__(data,
                                         names=columns,
                                         meta=meta)


    def __add__( self, other ):
        """ Overload the '+' operator.
        """
        raise NotImplementedError("I'm not yet sure how to do this")


    def __sub__(self, other):
        """ Overload the - operator """
        raise NotImplementedError("I'm not yet sure how to do this")


    def __mul__(self, other):
        """ Overload the * operator """
        raise NotImplementedError("I'm not yet sure how to do this")


    def __div__(self, other):
        """ Overload the / operator """
        raise NotImplementedError("I'm not yet sure how to do this")


    def __str__(self):
        """Prettier representation of object instanct"""

        return "Lightcurve Object"


    def concatenate(self, other):
        """ Concatenate two lightcurves

        The arrays are
        concatenated and re-sorted in order of the MJD array.

        """
        raise NotImplementedError("I'm not yet sure how to do this")

        out_obj = LightCurve()
        '''
        out_obj.gross = np.concatenate( [self.gross, other.gross] )
        out_obj.flux = np.concatenate( [self.flux, other.flux] )
        out_obj.background = np.concatenate( [self.background,
                                              other.background] )
        out_obj.mjd = np.concatenate( [self.mjd, other.mjd] )
        out_obj.times = np.concatenate( [self.times, other.times] )
        out_obj.bins = np.concatenate( [self.bins, other.bins] )

        sorted_index = np.argsort( out_obj.mjd )

        out_obj.gross = out_obj.gross[ sorted_index ]
        out_obj.flux = out_obj.flux[sorted_index]
        out_obj.background = out_obj.background[ sorted_index ]
        out_obj.mjd = out_obj.mjd[ sorted_index ]
        out_obj.times = out_obj.times[ sorted_index ]
        out_obj.bins = out_obj.bins[ sorted_index ]
        '''
        return out_obj

    @property
    def counts(self):
        """ Calculate counts array

        Counts are calculated as the gross - background


        Returns
        -------
        counts : np.ndarray
            Array containing the calculated counts

        """

        if not len(self['gross']):
            return self['gross'].copy()
        else:
            return self['gross'] - self['background']


    @property
    def error(self):
        """ Calculate error array

        Error is calculated as sqrt( gross + background )

        Returns
        -------
        error : np.ndarray
            Array containing the calculated error

        """

        if not len(self['gross']):
            return self['gross'].copy()
        else:
            return np.sqrt(self['gross'] + self['background'])


    @property
    def flux_error(self):
        """ Estimate the error in the flux array

        Flux error is currently estimated as having the same magnitude
        relative to flux as the error has to counts.  Very simple approximation
        for now

        Returns
        -------
        flux_error : np.ndarray
            Array containing the calculated error in flux

        """

        if not len(self['gross']):
            return self['gross'].copy()
        else:
            return (self['error'] / self['counts']) *  self['flux']


    @property
    def net(self):
        """ Calculate net array

        Net is calculated as counts / time

        Returns
        -------
        net : np.ndarray
            Array containing calculated net

        """

        if not len(self['counts']):
            return self['counts'].copy()
        else:
            return self['counts'] / self['bins'].astype(np.float64)


    @property
    def signal_to_noise(self):
        """ Quick signal to noise estimate

        S/N is calculated as the ratio of gross to error

        Returns
        -------
        sn_ratio : np.ndarray
            Array containing calculated signal to noise ratio

        """

        if not len(self['gross']):
            return self['gross'].copy()
        else:
            return self['gross'] / self['error']


    def filter(self, indices):
        """ Keep only a section of the whole lightcurve

        Modifies the LC in place

        Parameters
        ----------
        indices : np.ndarray
            Array of indexes to keep

        """

        self.gross = self.gross[indices]
        self.flux = self.flux[indices]
        self.background = self.background[indices]
        self.mjd = self.mjd[indices]
        self.times = self.times[indices]
        self.bins = self.bins[indices]


    def write(self, outname=None, clobber=False):
        """ Write lightcurve out to FITS file

        Parameters
        ----------
        outname : bool or str
            Either True/False, or output name

        clobber : bool
            Allow overwriting of existing file with same name

        """

        if isinstance(outname, str):
            self.outname = outname

        hdu_out = fits.HDUList(fits.PrimaryHDU())

        #-- Primary header
        hdu_out[0].header['GEN_DATE'] = (str(datetime.now()), 'Creation Date')
        hdu_out[0].header['LC_VER'] = (__version__, 'lightcurve version used')
        hdu_out[0].header['AP_VER'] = (astropy.__version__, 'Astropy version used')
        hdu_out[0].header['NP_VER'] = (np.__version__, 'Numpy version used')
        hdu_out[0].header['SP_VER'] = (scipy.__version__, 'Scipy version used')
        hdu_out[0].header.add_blank('', before='GEN_DATE')
        hdu_out[0].header.add_blank('{0:{fill}{align}67}'.format(' Lightcurve extraction keywords ',
                                                                 fill='-',
                                                                 align='^'), before='GEN_DATE')
        hdu_out[0].header.add_blank('', before='GEN_DATE')




        hdu_out[0].header.add_blank('', after='SP_VER')
        hdu_out[0].header.add_blank('{0:{fill}{align}67}'.format(' Source Data keywords ',
                                                                 fill='-',
                                                                 align='^'), after='SP_VER')
        hdu_out[0].header.add_blank('', after='SP_VER')

        try:
            hdu_out[0].header.extend(self.hdu[0].header, end=True)
        except AttributeError:
            pass


        #-- Ext 1 table data
        bins_col = fits.Column('bins', 'D', 'second', array=self.bins)
        times_col = fits.Column('times', 'D', 'second', array=self.times)
        mjd_col = fits.Column('mjd', 'D', 'MJD', array=self.mjd)
        gross_col = fits.Column('gross', 'D', 'counts', array=self.gross)
        counts_col = fits.Column('counts', 'D', 'counts', array=self.counts)
        net_col = fits.Column('net', 'D', 'counts/s', array=self.net)
        flux_col = fits.Column('flux', 'D', 'ergs/s', array=self.flux)
        flux_error_col = fits.Column('flux_error', 'D', 'ergs/s',
                                       array=self.flux_error)
        bkgnd_col = fits.Column('background', 'D', 'cnts',
                                  array=self.background)
        error_col = fits.Column('error', 'D', 'counts', array=self.error)

        tab = fits.new_table([bins_col,
                                times_col,
                                mjd_col,
                                gross_col,
                                counts_col,
                                net_col,
                                flux_col,
                                flux_error_col,
                                bkgnd_col,
                                error_col])
        hdu_out.append(tab)



        #-- Ext 1 header
        hdu_out[1].header.add_blank('')
        hdu_out[1].header.add_blank('{0:{fill}{align}67}'.format(' Lightcurve extraction keywords ',
                                                                 fill='-',
                                                                 align='^'))
        hdu_out[1].header.add_blank('')


        hdu_out[1].header.add_blank('')
        hdu_out[1].header.add_blank('{0:{fill}{align}67}'.format(' Source Data keywords ',
                                                                 fill='-',
                                                                 align='^'))
        hdu_out[1].header.add_blank('')

        try:
            hdu_out[1].header.extend(self.hdu[1].header, end=True)
        except AttributeError:
            pass


        if outname.endswith('.gz'):
            print("Nope, can't write to gzipped files")
            self.outname = self.outname[:-3]

        hdu_out.writeto(self.outname, clobber=clobber)

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

    print("Creating composite lightcurve from:")
    print("\n".join(filelist))

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
                other_file = [item for item in get_both_filenames(filename) if not item == filename][0]
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
            print(max_wave, min_wave)

            if max_wave < wmax:
                wmax = max_wave
            if min_wave > wmin:
                wmin = min_wave


    print('Using wavelength range of: {}-{}'.format(wmin, wmax))

    out_lc = lightcurve.LightCurve()

    for filename in filelist:
        new_lc = io.open(filename=filename,
                         step=1,
                         wlim=(wmin, wmax),
                         alt=None,
                         filter=True)

        out_lc = out_lc.concatenate(new_lc)

    out_lc.write(output, clobber=True)
