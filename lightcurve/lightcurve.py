"""
Selection of functions dealing with time-tag data in FITS format

"""

__all__ = ['LightCurve']

from astropy.utils.console import ProgressBar
from astropy.io import fits as pyfits
import os
import numpy as np

#-------------------------------------------------------------------------------

class LightCurve(object):
    """
    Turn an event list (*_rawtag_*.fits, *_corrtag_*.fits) into a lightcurve.
    
    Parameters
    ----------
    filename : string
        Name of fits file
    step : int
        Timestep for datapoints
    xlim : tuple, optional
        Lower and upper x bounds to extract
    ylim : tuple, optional
        Lower and upper y bounds to extract
    extract: bool, optional
        Extract from locations given in XTRACTAB instead of xlim,ylim
    fluxcal: bool, optional
        Convert total counts into flux using PHOTTAB
    normalize: bool, optional
        Normalize each lightcurve to mean value

    Returns
    -------
    times : array
        Array of time samples
    counts : array
        Array of counts or flux in each time sample
    error : array
        Array of error estimates for count array

    Examples
    --------
    >>> lightcurve( 'ipppssoot_corrtag.fits' )
    array([1, 2, 3, 4]), array([25, 25, 25, 25]), array([5, 5, 5, 5])

    """

    def __init__(self, **kwargs):
        """ Initialize and extract lightcurve from input corrtag
        """

        if 'filename' in kwargs:
            filetype = self._check_filetype( kwargs['filename'] )

            if filetype == 'corrtag':
                self.read_cos( kwargs['filename'] )
                self.extract()
            elif filetype == 'lightcurve':
                self.open_lightcurve( kwargs['filename'] )

        else:
            self.times = np.array( [] )
            self.mjd = np.array( [] )
            self.gross = np.array( [] )
            self.background = np.array( [] )


    def __add__( self, other ):
        """ Overload the '+' operator to concatenate lightcurve objects
        
        """

        out_obj = LightCurve()

        if isinstance( other, LightCurve ):
            out_obj.gross = np.concatenate( [self.gross, other.gross] )
            out_obj.background = np.concatenate( [self.background, other.background] )
            out_obj.mjd = np.concatenate( [self.mjd, other.mjd] )
            out_obj.times = np.concatenate( [self.times, other.times] )

            sorted_index = np.argsort( out_obj.mjd )

            out_obj.gross = out_obj.gross[ sorted_index ]
            out_obj.background = out_obj.background[ sorted_index ]
            out_obj.mjd = out_obj.mjd[ sorted_index ]
            out_obj.times = out_obj.times[ sorted_index ]

        else:
            out_obj.gross = self.gross + other
            out_obj.background = self.background + other
            out_obj.times = self.times
            out_obj.mjd = self.mjd
        
        return out_obj


    def __mul__( self, other ):
        """ Overload the * operator """
        
        out_obj = LightCurve()

        out_obj.gross = self.gross * other
        out_obj.background = self.background * other
        out_obj.times = self.times
        out_obj.mjd = self.mjd

        return out_obj


    def __div__( self, other ):
        """ Overload the / operator """
        
        out_obj = LightCurve()

        out_obj.gross = self.gross / other
        out_obj.background = self.background / other
        out_obj.times = self.times
        out_obj.mjd = self.mjd

        return out_obj


    def __str__(self):
        """Prettier representation of object instanct"""
        return "COS Lightcurve Object of %s" % ( ','.join( self.input_list ) ) 
 
    @property
    def counts(self):
        """ Calculate counts array """
        
        if not len( self.gross ):
            return self.gross.copy()
        else:
            return self.gross - self.background


    @property
    def error(self):
        """ Calculate error array """

        if not len(self.counts):
            return self.counts.copy()
        else:
            return np.sqrt( self.counts + self.background )

        
    @property
    def net(self):
        """ Calculate net array """

        if not len(self.counts):
            return self.counts.copy()
        else:
            return self.counts / self.times.astype( np.float64 )


    @classmethod
    def _check_filetype(self, filename):
        """Determine the type of data being input"""
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

        lightcurve_names = set( ['TIME',
                                 'MJD',
                                 'COUNTS',
                                 'NET',
                                 'FLUX',
                                 'BACKGROUND',
                                 'ERROR'] )

        hdu = pyfits.open( filename )
        input_names = set( [item.upper() for 
                           item in hdu[1].data.names ] )
        hdu.close()

        if input_names == corrtag_names:
            filetype = 'corrtag'
        elif input_names == lightcurve_list:
            filetype = 'lightcurve'
        else:
            filetype = None

        return filetype


    def read_cos(self, filename ):
        """ Extract light curve from COS data"""

        self.hdu = pyfits.open( filename )
        self.input_filename = filename

        self._get_hdus()


    @classmethod
    def open_lightcurve(cls, filename):    
        """ Read fits lightcurve from fits file back into object"""

        out_obj = cls()
        
        hdu = pyfits.open( filename)
        
        out_obj.times = hdu[1].data['times']
        out_obj.gross = hdu[1].data['gross']
        out_obj.mjd = hdu[1].data['mjd']
        #out_obj.flux = hdu[1].data['flux']
        out_obj.background = hdu[1].data['background']

        return out_obj


    def normalize(self):
        """ Normalize arrays around mean"""
        self.gross = self.gross / self.gross.mean()
        self.background = self.background / self.background.mean()


    def extract(self, step=1, xlim=(0,16384), wlim=(-1,10000), ylim=(0,1024) ):
        """ Loop over HDUs and extract the lightcurve

        """
        SECOND_PER_MJD = 1.15741e-5
        end = min( self.hdu['events'].data[ 'time' ].max(), self.hdu[1].header['EXPTIME'] )
        all_steps = np.arange(0, end, step)[:-1]

        gross = 0
        background = 0
        for segment, hdu in self.hdu_dict.iteritems():
            ystart, yend = self._get_extraction_region( hdu, segment, 'spectrum' )
            index = self._extract_index(hdu,
                                        xlim[0], xlim[1], 
                                        ystart, yend, 
                                        wlim[0], wlim[1],
                                        hdu[1].header['sdqflags'] )
            gross += np.histogram( hdu[ 'events' ].data['time'][index], all_steps )[0]

            bstart, bend = self._get_extraction_region( hdu, segment, 'background1' )
            index = self._extract_index(hdu,
                                         xlim[0], xlim[1], 
                                         bstart, bend, 
                                         wlim[0], wlim[1],
                                         hdu[1].header['sdqflags'] )
            b_gross = np.histogram( hdu[ 'events' ].data['time'][index], all_steps )[0]

            b_correction = ( (bend - bstart) / (yend -ystart) )
            background += b_gross * b_correction

        self.gross = gross
        self.background = background
        self.mjd = self.hdu[1].header[ 'EXPSTART' ] + np.array( all_steps[:-1] ) * SECOND_PER_MJD
        self.times = np.ones( len( gross ) ) * step


    def _check_output(self, writeto):
        """ Determine what was supplied to writeto and set self.outname

        If a string, lightcurve will be written to value.
        If True, output will be rootname_curve.fits

        """

        if writeto:
            if isinstance( writeto, bool ):
                self.outname = self.hdu[0].header['rootname'] + '_curve.fits'
            elif isinstance( writeto, str ):
                self.outname = writeto
            else:
                raise NameError( 'writeto should be True/False our string :' )
            
            if os.path.exists( self.outname ) and not self.clobber:
                raise IOError( 'Output already exists %s'% self.outname )
        else:
            pass


    def _get_hdus(self):
        """Populate HDU dictionary from available corrtag files

        """

        if self.hdu[0].header['detector'] == 'NUV':
            all_filenames = [ self.input_filename ]

            if self.hdu[0].header['opt_elem'] == 'G230L':
                hdu_dict = {'A':self.hdu, 'B':self.hdu}
            else:
                hdu_dict = {'A':self.hdu, 'B':self.hdu, 'C':self.hdu}

        elif self.hdu[0].header['detector'] == 'FUV':
            file_a, file_b = self._get_both_filenames( self.input_filename )

            hdu_dict = dict()
            all_filenames = []
            for name, letter in zip( [file_a, file_b], ['A', 'B'] ):
                if os.path.exists( name ):
                    all_filenames.append( name )
                    hdu_dict[ letter ] = pyfits.open( name, ignore_missing_end=True )

        self.hdu_dict = hdu_dict
        self.input_list = all_filenames


    def _get_both_filenames(self, filename ):
        """ Get a list of both filenames for FUV data

        Regardless if rootname_corrtag_a.fits or rootname_corrtag_b.fits 
        is passed in, both will be returned in a list.

        """

        assert pyfits.getval( filename, 'DETECTOR' ) == 'FUV', 'This only works for FUV data'

        if pyfits.getval(filename, 'SEGMENT') == 'FUVA':
            other_filename = filename.replace('_a.fits', '_b.fits')
        elif pyfits.getval(filename, 'SEGMENT') == 'FUVB':
            other_filename = filename.replace('_b.fits', '_a.fits')

        filename_list = [filename, other_filename]
        filename_list.sort()

        return ( filename_list[0], filename_list[1] )


    def _get_extraction_region(self, hdu, segment, mode='spectrum' ):
        """Get y_start,y_end for given extraction

        """

        if mode == 'spectrum':
            height =  hdu[1].header['SP_HGT_%s' % (segment)]
            location = hdu[1].header['SP_LOC_%s' % (segment)]

        elif mode == 'background1':
            height =  hdu[1].header['B_HGT1_%s' % (segment)]
            location = hdu[1].header['B_BKG1_%s' % (segment)] 

        elif mode == 'background2':
            height =  hdu[1].header['B_HGT2_%s' % (segment)]
            location = hdu[1].header['B_BKG2_%s' % (segment)] 

        y_start = location - height/2
        y_end = location + height/2

        return y_start, y_end


    def _extract_index(self, hdu, x_start, x_end, 
                        y_start, y_end, w_start, w_end, sdqflags=0):
        """
        Extract counts from given HDU using input parameters.

        Lyman Alpha line counts will be excluded by default.

        """

        data_index = np.where( ( hdu[1].data['XCORR'] >= x_start ) & 
                               ( hdu[1].data['XCORR'] < x_end ) &

                               ( hdu[1].data['YCORR'] >= y_start ) & 
                               ( hdu[1].data['YCORR'] < y_end ) &

                               ~( hdu[1].data['DQ'] & sdqflags ) &

                               ( (hdu[1].data['WAVELENGTH'] > w_start) & 
                                 (hdu[1].data['WAVELENGTH'] < w_end) ) &

                               ( (hdu[1].data['WAVELENGTH'] > 1217) | 
                                 (hdu[1].data['WAVELENGTH'] < 1214) ) &
                               ( (hdu[1].data['WAVELENGTH'] > 1308) | 
                                 (hdu[1].data['WAVELENGTH'] < 1300) ) )[0]

        return data_index


    def _get_flux_correction(self, fluxtab, opt_elem, cenwave, 
                             aperture, minwave, maxwave):
        """ Return integrated flux calibration over given wavelength range 
        """

        if '$' in fluxtab:
            fluxpath, fluxfile = fluxtab.split( '$' )
            try:
                fluxfile = os.path.join( os.environ['lref'], fluxfile )
            except KeyError:
                fluxfile = None
        else:
            fluxpath = './'
            fluxfile = fluxtab

        if (not fluxfile) or (not os.path.exists( fluxfile ) ):
            print ' WARNING: Fluxfile not available %s,' % fluxfile
            print ' using unity flux calibration instead.'
            return 1

        flux_hdu = pyfits.open( fluxfile )
        flux_index = np.where( (flux_hdu[1].data['OPT_ELEM'] == opt_elem) &
                               (flux_hdu[1].data['CENWAVE'] == cenwave) &
                               (flux_hdu[1].data['APERTURE'] == aperture) )[0]

        mean_response_list = []
        for sens_curve, wave_curve in zip( flux_hdu[1].data[flux_index]['SENSITIVITY'],
                                          flux_hdu[1].data[flux_index]['WAVELENGTH'] ):
            wave_index = np.where( (wave_curve >= minwave) &
                                   (wave_curve <= maxwave ) )[0]

            if not len(wave_index): continue

            mean_response_list.append( sens_curve[ wave_index ].sum() )

        total_fluxcorr = np.sum( mean_response_list )

        return total_fluxcorr


    def write(self, outname=None, clobber=False):
        """ Write out to FITS file
        """

        if isinstance( outname, str ):
            self.outname = outname

        hdu_out = pyfits.HDUList(pyfits.PrimaryHDU())

        try: hdu_out[0].header = self.hdu[0].header
        except: pass 

        times_col = pyfits.Column('times', 'D', 'second', array=self.times)
        mjd_col = pyfits.Column('mjd', 'D', 'MJD', array=self.mjd) 
        gross_col = pyfits.Column('gross', 'D', 'counts', array=self.counts)    
        counts_col = pyfits.Column('counts', 'D', 'counts', array=self.counts)
        net_col = pyfits.Column('net', 'D', 'counts/s', array=self.net)
        #flux_col = pyfits.Column('flux', 'D', 'ergs/s', array=self.flux)
        bkgnd_col = pyfits.Column('background', 'D', 'cnts', array=self.background)
        error_col = pyfits.Column('error', 'D', 'counts', array=self.error)
        
        tab = pyfits.new_table( [time_col,
                                 mjd_col,
                                 gross_col,
                                 counts_col,
                                 net_col,
                                 bkgnd_col,
                                 error_col] )
        hdu_out.append( tab )

        hdu_out.writeto( self.outname, clobber=clobber)  

