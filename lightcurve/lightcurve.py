"""
Class and tools for extracting lightcurve data from time-tag event lists

"""

__all__ = ['LightCurve']

from astropy.io import fits as pyfits
import numpy as np

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

    def __init__(self):
        """ Initialize and extract lightcurve from input corrtag

        lightcurves must contain:
            gross
            times (timestep)
            mjd
            background
            flux

        """

        self.gross = np.array( [] )
        self.times = np.array( [] )
        self.mjd = np.array( [] )
        self.background = np.array( [] )
        self.flux = np.array( [] )
        self.times = np.array( [] )


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

        if isinstance( other, LightCurve ):
            out_obj.gross = self.gross / other.gross
            out_obj.flux = self.flux / other.gross
            out_obj.background = self.background / other.background
            out_obj.times = self.times / other.times
            out_obj.mjd = self.mjd

        else:
            out_obj.gross = self.gross / other
            out_obj.flux = self.flux / other
            out_obj.background = self.background / other
            out_obj.times = self.times
            out_obj.mjd = self.mjd

        return out_obj


    def __str__(self):
        """Prettier representation of object instanct"""
        
        return "Lightcurve Object"
 

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

        if not len(self.gross):
            return self.gross.copy()
        else:
            return np.sqrt( self.gross + self.background )


    @property
    def flux_error(self):
        """ Estimate the error in the flux array 
        
        Flux error is currently estimated as having the same magnitude
        relative to flux as the error has to counts.  Very simple approximation
        for now

        """

        if not len(self.gross):
            return self.gross.copy()
        else:
            return (self.error / self.counts ) *  self.flux


    @property
    def net(self):
        """ Calculate net array """

        if not len(self.counts):
            return self.counts.copy()
        else:
            return self.counts / self.times.astype( np.float64 )

        
    @property
    def signal_to_noise(self):
        """ Quick signal to noise estimate
        """

        if not len(self.gross):
            return self.gross.copy()
        else:
            return self.gross / self.error


    def normalize(self):
        """ Normalize arrays around mean"""
        self.gross = self.gross / self.gross.mean()
        self.background = self.background / self.background.mean()


    def write(self, outname=None, clobber=False):
        """ Write out to FITS file
        """

        if isinstance( outname, str ):
            self.outname = outname

        hdu_out = pyfits.HDUList(pyfits.PrimaryHDU())

        try: hdu_out[0].header = self.hdu[0].header
        except AttributeError: pass 

        times_col = pyfits.Column('times', 'D', 'second', array=self.times)
        mjd_col = pyfits.Column('mjd', 'D', 'MJD', array=self.mjd) 
        gross_col = pyfits.Column('gross', 'D', 'counts', array=self.gross)    
        counts_col = pyfits.Column('counts', 'D', 'counts', array=self.counts)
        net_col = pyfits.Column('net', 'D', 'counts/s', array=self.net)
        flux_col = pyfits.Column('flux', 'D', 'ergs/s', array=self.flux)
        flux_error_col = pyfits.Column('flux_error', 'D', 'ergs/s', 
                                       array=self.flux_error)
        bkgnd_col = pyfits.Column('background', 'D', 'cnts', 
                                  array=self.background)
        error_col = pyfits.Column('error', 'D', 'counts', array=self.error)
        
        tab = pyfits.new_table( [times_col,
                                 mjd_col,
                                 gross_col,
                                 counts_col,
                                 net_col,
                                 flux_col,
                                 flux_error_col,
                                 bkgnd_col,
                                 error_col] )
        hdu_out.append( tab )

        if self.outname.endswith('.gz'):
            print "Nope, can't write to gzipped files"
            self.outname = self.outname[:-3]

        hdu_out.writeto( self.outname, clobber=clobber)  

