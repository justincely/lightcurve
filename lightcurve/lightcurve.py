"""
Selection of functions dealing with time-tag data in FITS format

"""

__all__ = ['LightCurve']

from astropy.utils.console import ProgressBar
from astropy.io import fits as pyfits
import os
import numpy as np
from scipy.interpolate import interp1d

#-------------------------------------------------------------------------------

class LightCurve(object):
    """
    Turn an event list (*_rawtag_*.fits, *_corrtag_*.fits) into a lightcurve.
    
    Parameters
    ----------
    filename : string
        Name of fits file

    Returns
    -------
    lightcurve object

    Examples
    --------
    >>> obj = LightCurve( filename='ipppssoot_corrtag.fits' )
    obj.gross
    obj.background
    obj.error

    """

    def __init__(self, **kwargs):
        """ Initialize and extract lightcurve from input corrtag
        """

        if 'step' in kwargs: step=kwargs['step']
        else: step = 1

        if 'wlim' in kwargs: wlim = kwargs['wlim']
        else: wlim = (1,10000)

        if 'filename' in kwargs:
            filetype = self._check_filetype( kwargs['filename'])

            if filetype == 'corrtag':
                self.read_cos( kwargs['filename'])
                self.extract( step=step, wlim=wlim )
            elif filetype == 'lightcurve':
                self.open_lightcurve( kwargs['filename'] )

        else:
            self.times = np.array( [] )
            self.flux = np.array( [] )
            self.mjd = np.array( [] )
            self.gross = np.array( [] )
            self.background = np.array( [] )


    def __add__( self, other ):
        """ Overload the '+' operator.  If other is a float or int, the entire
        the value is added to the gross and background arrays.  If other is
        another LightCurve object, the arrays are concatenated and re-sorted
        in order of the MJD array.
        
        """

        out_obj = LightCurve()

        if isinstance( other, LightCurve ):
            out_obj.gross = np.concatenate( [self.gross, other.gross] )
            out_obj.flux = np.concatenate( [self.flux, other.flux] )
            out_obj.background = np.concatenate( [self.background, 
                                                  other.background] )
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
        out_obj.flux = self.flux * other
        out_obj.background = self.background * other
        out_obj.times = self.times
        out_obj.mjd = self.mjd

        return out_obj


    def __div__( self, other ):
        """ Overload the / operator """
        
        out_obj = LightCurve()

        out_obj.gross = self.gross / other
        out_obj.flux = self.flux / other
        out_obj.background = self.background / other
        out_obj.times = self.times
        out_obj.mjd = self.mjd

        return out_obj


    def __str__(self):
        """Prettier representation of object instanct"""
        
        if not 'input_list' in dir( self ):
            self.input_list = ['None']

        return "Lightcurve Object of %s" % ( ','.join( self.input_list ) ) 
 

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
        """Determine the type of data being input.
        
        File type is determined by the culumns in the first data extension.

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

        lightcurve_names = set( ['TIME',
                                 'MJD',
                                 'GROSS',
                                 'COUNTS',
                                 'NET',
                                 'BACKGROUND',
                                 'ERROR'] )

        hdu = pyfits.open( filename )
        input_names = set( [item.upper() for 
                           item in hdu[1].data.names ] )
        hdu.close()

        if input_names == corrtag_names:
            filetype = 'corrtag'
        elif input_names == lightcurve_names:
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
        """ Read lightcurve from fits file back into object"""

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


    def extract(self, step=1, xlim=(0,16384), wlim=(2,10000), ylim=None ):
        """ Loop over HDUs and extract the lightcurve

        """


        SECOND_PER_MJD = 1.15741e-5
        end = min( self.hdu['events'].data[ 'time' ].max(), self.hdu[1].header['EXPTIME'] )
        all_steps = np.arange(0, end+step, step)

        if all_steps[-1] > end:
            truncate = True
        else:
            truncate = False

        gross = 0
        flux = 0
        background = 0
        background_flux = 0

        for segment, hdu in self.hdu_dict.iteritems():

            ### Source Calculation
            
            if not ylim:
                ystart, yend = self._get_extraction_region( hdu, segment, 'spectrum' )
            else:
                ystart, yend = ylim[0], ylim[1]

            index = self._extract_index(hdu,
                                        xlim[0], xlim[1], 
                                        ystart, yend, 
                                        wlim[0], wlim[1],
                                        hdu[1].header['sdqflags'] )
            gross += np.histogram( hdu[ 'events' ].data['time'][index], all_steps, 
                                   weights=hdu[ 'events' ].data['epsilon'][index] )[0]

            response_array = self._get_fluxes( hdu, index )
            tds_corr = self._get_tds( hdu, index )

            weights = hdu[ 'events' ].data['epsilon'][index] / step / tds_corr / response_array 
            flux +=  np.histogram( hdu[ 'events' ].data['time'][index], all_steps, weights=weights)[0]


            ### Background calculation

            bstart, bend = self._get_extraction_region( hdu, segment, 'background1' )
            index = self._extract_index(hdu,
                                        xlim[0], xlim[1], 
                                        bstart, bend, 
                                        wlim[0], wlim[1],
                                        hdu[1].header['sdqflags'] )
            
            b_corr = ( (bend - bstart) / (yend -ystart) )
            background += b_corr * np.histogram( hdu[ 'events' ].data['time'][index], all_steps, 
                                                 weights=hdu[ 'events' ].data['epsilon'][index]  )[0]
            
            response_array = self._get_fluxes( hdu, index )
            tds_corr = self._get_tds( hdu, index )
            weights = hdu[ 'events' ].data['epsilon'][index] / step / tds_corr / response_array
            background_flux +=  b_corr * np.histogram( hdu[ 'events' ].data['time'][index], all_steps, weights=weights )[0] 


        self.gross = gross
        self.flux = flux - background_flux
        self.background = background
        self.mjd = self.hdu[1].header[ 'EXPSTART' ] + np.array( all_steps[:-1] ) * SECOND_PER_MJD
        self.times = np.ones( len( gross ) ) * step

        if truncate:
            self.gross = self.gross[:-1]
            self.flux = self.flux[:-1]
            self.background = self.background[:-1]
            self.mjd = self.mjd[:-1]
            self.times = self.times[:-1]


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

        #if hdu[0].header['OPT_ELEM'] == 'G140L':
        #    w_start = max( 920, wstart )
        #    w_end = min( 1800, wend )

        data_index = np.where( ( hdu[1].data['XCORR'] >= x_start ) & 
                               ( hdu[1].data['XCORR'] < x_end ) &

                               ( hdu[1].data['YCORR'] >= y_start ) & 
                               ( hdu[1].data['YCORR'] < y_end ) &

                               np.logical_not( hdu[1].data['DQ'] & sdqflags ) &

                               ( (hdu[1].data['WAVELENGTH'] > w_start) & 
                                 (hdu[1].data['WAVELENGTH'] < w_end) ) &

                               ( (hdu[1].data['WAVELENGTH'] > 1217) | 
                                 (hdu[1].data['WAVELENGTH'] < 1214) ) &
                               ( (hdu[1].data['WAVELENGTH'] > 1308) | 
                                 (hdu[1].data['WAVELENGTH'] < 1300) ) )[0]

        return data_index


    def _get_tds(self, hdu, index):
        """ return tds correction for each event
        """

        try:
            tdstab = hdu[0].header['TDSTAB']
        except KeyError:
            tdstab = ''


        if '$' in tdstab:
            tdspath, tdsfile = tdstab.split( '$' )
            try:
                tdsfile = os.path.join( os.environ[tdspath], tdsfile )
            except KeyError:
                tdsfile = None
        else:
            tdspath = './'
            tdsfile = tdstab


        if (not tdsfile) or (not os.path.exists( tdsfile ) ):
            print ' WARNING: Fluxfile not available %s,' % tdsfile
            print ' using unity flux calibration instead.'
            return np.ones( hdu[ 'events' ].data['time'].shape )[index]
        
        tds_data = pyfits.getdata( tdsfile,ext=1 )
        REF_TIME = pyfits.getval( tdstab,'REF_TIME',ext=1)
        mode_index = np.where( (tds_data['OPT_ELEM'] == hdu[0].header['opt_elem'] ) &
                               (tds_data['APERTURE'] == hdu[0].header['aperture'] ) &
                               (tds_data['SEGMENT'] == hdu[0].header['segment'] ) )[0]

        mode_line = tds_data[ mode_index ]

        tds_nwl = mode_line['NWL'][0]
        tds_nt = mode_line['NT'][0]
        tds_wavelength = mode_line['WAVELENGTH'][0]
        smaller_index = np.where( mode_line['TIME'][0][:tds_nt] < hdu[1].header['EXPSTART'] )[0]
        time_index = smaller_index.max()

        tds_slope = mode_line['SLOPE'][0]
        tds_intercept = mode_line['INTERCEPT'][0]

        calculate_drop = lambda time, ref_time, slope, intercept: \
            ( (( time - ref_time ) * slope) / (365.25*100) ) + intercept

        correction_array = np.zeros( len(tds_wavelength) )
        for i,wave in enumerate( tds_wavelength ):
            correction = calculate_drop( hdu[1].header['EXPSTART'], 
                                         REF_TIME, 
                                         tds_slope[time_index,i], 
                                         tds_intercept[time_index,i] )
            correction_array[i] = correction

        #-------
        # interpolate onto input arrays
        #------

        interp_function = interp1d( tds_wavelength, correction_array, 1 )
        response_correction = interp_function( hdu[1].data['wavelength'][index] )

        return response_correction


    def _get_fluxes(self, hdu, index):
        """ Return response curve and wavelength range for specified mode
        """

        try:
            fluxtab = hdu[0].header['FLUXTAB']
        except KeyError:
            fluxtab = ''


        if '$' in fluxtab:
            fluxpath, fluxfile = fluxtab.split( '$' )
            try:
                fluxfile = os.path.join( os.environ[fluxpath], fluxfile )
            except KeyError:
                fluxfile = None
        else:
            fluxpath = './'
            fluxfile = fluxtab

        if (not fluxfile) or (not os.path.exists( fluxfile ) ):
            print ' WARNING: Fluxfile not available %s,' % fluxfile
            print ' using unity flux calibration instead.'
            return np.ones( hdu[ 'events' ].data['time'].shape )[index]


        flux_hdu = pyfits.open( fluxfile )
        flux_index = np.where( (flux_hdu[1].data['SEGMENT'] == hdu[0].header['segment']) &
                               (flux_hdu[1].data['OPT_ELEM'] == hdu[0].header['opt_elem']) &
                               (flux_hdu[1].data['CENWAVE'] == hdu[0].header['cenwave']) &
                               (flux_hdu[1].data['APERTURE'] == hdu[0].header['aperture']) )[0]


        interp_func = interp1d( flux_hdu[1].data[flux_index]['WAVELENGTH'].flatten(), 
                                flux_hdu[1].data[flux_index]['SENSITIVITY'].flatten() )

        all_resp = interp_func( hdu[ 'events' ].data['wavelength'][index] )

        return all_resp


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
        gross_col = pyfits.Column('gross', 'D', 'counts', array=self.gross)    
        counts_col = pyfits.Column('counts', 'D', 'counts', array=self.counts)
        net_col = pyfits.Column('net', 'D', 'counts/s', array=self.net)
        flux_col = pyfits.Column('flux', 'D', 'ergs/s', array=self.flux)
        bkgnd_col = pyfits.Column('background', 'D', 'cnts', array=self.background)
        error_col = pyfits.Column('error', 'D', 'counts', array=self.error)
        
        tab = pyfits.new_table( [times_col,
                                 mjd_col,
                                 gross_col,
                                 counts_col,
                                 net_col,
                                 flux_col,
                                 bkgnd_col,
                                 error_col] )
        hdu_out.append( tab )

        hdu_out.writeto( self.outname, clobber=clobber)  


#--------------------------------------------------------------

def backout_tds( response, wavelength, mjd, opt_elem, aperture, segment, tdstab='/grp/hst/cdbs/lref/w7h1935dl_tds.fits'):
    """
    Backs out TDS from response curve

    #COPY from COS reference files repository
    """

    #print 'Removing TDS for:', opt_elem, aperture, segment
    #print 'Using:',tdstab,' From:',mjd


    #-------#
    # get correction array 
    #-------#
    tds_data = pyfits.getdata( tdstab,ext=1 )
    REF_TIME = pyfits.getval( tdstab,'REF_TIME',ext=1)
    mode_index = np.where( (tds_data['OPT_ELEM'] == opt_elem) &
                           (tds_data['APERTURE'] == aperture) &
                           (tds_data['SEGMENT'] == segment) )[0]

    mode_line = tds_data[ mode_index ]
   
    tds_nwl = mode_line['NWL'][0]
    tds_nt = mode_line['NT'][0]
    tds_wavelength = mode_line['WAVELENGTH'][0]
    smaller_index = np.where( mode_line['TIME'][0][:tds_nt] < mjd )[0]
    time_index = smaller_index.max()

    tds_slope = mode_line['SLOPE'][0]
    tds_intercept = mode_line['INTERCEPT'][0]

    correction_array = np.zeros( len(tds_wavelength) )
    for i,wave in enumerate( tds_wavelength ):
        correction = calculate_drop( mjd, REF_TIME, tds_slope[time_index,i], tds_intercept[time_index,i] )
        correction_array[i] = correction

    #-------
    # interpolate onto input arrays
    #------

    interp_function = interp1d( tds_wavelength, correction_array, 1 )
    interp_correction = interp_function( wavelength )

    return response / interp_correction

 #--------------------------------------------------------------

