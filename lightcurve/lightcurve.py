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

    def __init__(self):
        """ Initialize and extract lightcurve from input corrtag
        """
        self.times = np.array( [] )
        self.mjd = np.array( [] )
        self.counts = np.array( [] )
        self.background = np.array( [] )

        #self.net = np.array( [] )
        self.flux = np.array( [] )
        #self.error = np.array( [] )


    def __add__( self, other ):
        """ Overload the '+' operator to concatenate lightcurve objects
        
        """

        out_obj = LightCurve()

        out_obj.counts = np.concatenate( [self.counts, other.counts] )
        out_obj.flux = np.concatenate( [self.flux, other.flux] )
        out_obj.background = np.concatenate( [self.background, other.background] )
        out_obj.mjd = np.concatenate( [self.mjd, other.mjd] )
        out_obj.times = np.concatenate( [self.times, other.times] )

        sorted_index = np.argsort( out_obj.mjd )
        
        out_obj.counts = out_obj.counts[ sorted_index ]
        out_obj.flux = out_obj.flux[ sorted_index ]
        out_obj.background = out_obj.background[ sorted_index ]
        out_obj.mjd = out_obj.mjd[ sorted_index ]
        out_obj.times = out_obj.times[ sorted_index ]
        
        return out_obj


    def __mul__( self, other ):
        """ Overload the * operator """
        
        out_obj = LightCurve()

        out_obj.counts = self.counts * other
        out_obj.flux = self.flux * other
        out_obj.background = self.background * other
        out_obj.times = self.times
        out_obj.mjd = self.mjd

        return out_obj


    def __div__( self, other ):
        """ Overload the / operator """
        
        out_obj = LightCurve()

        out_obj.counts = self.counts / other
        out_obj.flux = self.flux / other
        out_obj.background = self.background / other
        out_obj.times = self.times
        out_obj.mjd = self.mjd

        return out_obj


    def __str__(self):
        """Prettier representation of object instanct"""
        return "COS Lightcurve Object of %s" % ( ','.join( self.input_list ) ) 
 

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
    def check_filetype(self, filename):
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

    @classmethod
    def extract_from_cos(cls, filename, step=5, xlim=None, wlim=(-1, 10000), 
                 fluxtab=None, normalize=False, writeto=None, clobber=False):
        """ Extract light curve from COS data"""

        out_obj = cls()

        out_obj.hdu = pyfits.open( filename )
        out_obj.input_filename = filename
        out_obj.clobber = clobber

        out_obj._check_output( writeto )

        out_obj.fluxtab = fluxtab or out_obj.hdu[0].header[ 'FLUXTAB' ]
        out_obj.detector = out_obj.hdu[0].header[ 'DETECTOR' ]
        out_obj.opt_elem = out_obj.hdu[0].header[ 'OPT_ELEM' ]
        out_obj.cenwave = out_obj.hdu[0].header[ 'CENWAVE' ]
        out_obj.aperture = out_obj.hdu[0].header[ 'APERTURE' ]
        out_obj.sdqflags = out_obj.hdu[1].header[ 'SDQFLAGS' ]

        out_obj._get_hdus()

        print 'Making lightcurve from: ' + ', '.join( out_obj.input_list )
        print 'Extracting on stripes: ' +','.join( np.sort( out_obj.hdu_dict.keys()) )
        print 'DETECTOR: %s'% out_obj.detector
        print 'APERTURE: %s'% out_obj.aperture
        print 'OPT_ELEM: %s'% out_obj.opt_elem
        print 'CENWAVE : %s'% out_obj.cenwave

        if (not xlim) and (out_obj.detector == 'FUV'):
            xlim = (0, 16384)
        elif (not xlim) and (out_obj.detector == 'NUV'):
            xlim = (0, 1024)

        out_obj.generate_lightcurve(xlim, wlim, step)

        if normalize: out_obj.normalize()

        if writeto:
            out_obj.write( clobber=out_obj.clobber  )

        return out_obj

    @classmethod
    def open_lightcurve(cls, filename):    
        """ Read fits lightcurve from fits file back into object"""
        out_obj = cls()
        
        hdu = pyfits.open( filename)
        
        out_obj.times = hdu[1].data['time']
        out_obj.mjd = hdu[1].data['mjd']
        out_obj.counts = hdu[1].data['counts']
        out_obj.net = hdu[1].data['net']
        out_obj.flux = hdu[1].data['flux']
        out_obj.background = hdu[1].data['background']
        out_obj.error = hdu[1].data['error']

        return out_obj

    def normalize(self):
        """ Normalize arrays around mean"""
        self.error = self.error / self.counts
        self.counts = self.counts / self.counts.mean()
        self.background = self.background / self.background.mean()
        self.flux = self.flux / self.flux.mean()


    def generate_lightcurve(self, xlim, wlim, step):
        """ Loop over HDUs and extract the lightcurve
        
        """
        
        SECOND = 1.15741e-5
        EXPSTART = self.hdu[1].header[ 'EXPSTART' ]
        end = min( self.hdu['events'].data[ 'time' ].max(), self.hdu[1].header['EXPEND'] )
        

        print "Total Time = %ds" %  ( int(end) )
        print "Bins = %ds" % (step)

        all_counts = []
        all_bkgnd = []
        all_steps = range(0, end, step)[:-1]
        all_mjd = []

        int_flux_curve = self._get_flux_correction( self.fluxtab, self.opt_elem, 
                                                    self.cenwave, self.aperture, 
                                                    wlim[0], wlim[1] )

        for start in ProgressBar.iterate( all_steps ):
            sub_count = []
            sub_bkgnd = []
            for segment, hdu in self.hdu_dict.iteritems():
                ystart, yend = self._get_extraction_region( hdu, segment, 'spectrum' )
                counts, wmin, wmax = self._extract_counts(hdu, start, 
                                                          start + step, 
                                                          xlim[0], xlim[1], 
                                                          ystart, yend, 
                                                          wlim[0], wlim[1],
                                                          self.sdqflags )

                bstart, bend = self._get_extraction_region( hdu, segment, 'background1' )
                b_counts, b_wmin, b_wmax = self._extract_counts(hdu, start, 
                                                                start + step, 
                                                                xlim[0], xlim[1], 
                                                                bstart, bend, 
                                                                wlim[0], wlim[1],
                                                                self.sdqflags )

                b_correction = ( (bend - bstart) / (yend -ystart) )

                background = b_counts * b_correction

                sub_count.append( counts )
                sub_bkgnd.append( background )

            sample_counts = np.sum( sub_count )
            sample_bkgnd = np.sum( sub_bkgnd )

            all_counts.append( sample_counts - sample_bkgnd )
            all_bkgnd.append( sample_bkgnd )
            all_mjd.append( EXPSTART + (start * SECOND) )

        self.counts = np.array( all_counts )
        self.flux = self.net / int_flux_curve
        self.background = np.array( all_bkgnd )
        self.mjd = np.array( all_mjd )
        self.times = np.ones( len(self.counts) ) * step


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

        if self.detector == 'NUV':
            all_filenames = [ self.input_filename ]

            if self.opt_elem == 'G230L':
                hdu_dict = {'A':self.hdu, 'B':self.hdu}
            else:
                hdu_dict = {'A':self.hdu, 'B':self.hdu, 'C':self.hdu}

        elif self.detector == 'FUV':
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


    def _extract_counts(self, hdu, start, end, x_start, x_end, 
                        y_start, y_end, w_start, w_end, sdqflags=0):
        """
        Extract counts from given HDU using input parameters.

        Lyman Alpha line counts will be excluded by default.

        """

        data_index = np.where( ( hdu[1].data['TIME'] >= start ) & 
                               ( hdu[1].data['TIME'] < end ) &

                               ( hdu[1].data['XCORR'] >= x_start ) & 
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

        if len(data_index):
            minwave = hdu[1].data[ data_index ]['wavelength'].min()
            maxwave = hdu[1].data[ data_index ]['wavelength'].max()
        else:
            minwave = 0
            maxwave = 0

        return len(data_index), minwave, maxwave


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

        time_col = pyfits.Column('time', 'D', 'second', array=self.times)
        mjd_col = pyfits.Column('mjd', 'D', 'MJD', array=self.mjd)     
        counts_col = pyfits.Column('counts', 'D', 'counts', array=self.counts)
        net_col = pyfits.Column('net', 'D', 'counts/s', array=self.net)
        flux_col = pyfits.Column('flux', 'D', 'ergs/s', array=self.flux)
        bkgnd_col = pyfits.Column('background', 'D', 'cnts', array=self.background)
        error_col = pyfits.Column('error', 'D', 'counts', array=self.error)
        
        tab = pyfits.new_table( [time_col, mjd_col, counts_col, net_col, 
                                 flux_col, bkgnd_col, error_col] )
        hdu_out.append( tab )

        hdu_out.writeto( self.outname, clobber=clobber, output_verify='fix' )  

