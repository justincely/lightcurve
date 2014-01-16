"""
Contains the basic LightCurve object which is inherited by the subclasses.

"""

from __future__ import print_function

__all__ = ['LightCurve']

import astropy
import scipy
import numpy as np
from datetime import datetime

from .version import version as  __version__
from astropy.io import fits as pyfits


#-------------------------------------------------------------------------------

class LightCurve(object):
    """
    Returns
    -------
    lightcurve object

    """

    def __init__(self):
        """ Instantiate an empty LightCurve object.

        LightCurves must contain:
            gross
            times (timestep)
            mjd
            background
            flux

        All values are currently set to empty arrays.

        """

        self.gross = np.array( [] )
        self.times = np.array( [] )
        self.mjd = np.array( [] )
        self.background = np.array( [] )
        self.flux = np.array( [] )
        self.times = np.array( [] )


    def __add__( self, other ):
        """ Overload the '+' operator. 

        If other is another LightCurve object, the arrays are 
        concatenated and re-sorted in order of the MJD array.
        
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
            raise NotImplementedError("I'm not yet sure how to do this")
        
        return out_obj

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
 

    @property
    def counts(self):
        """ Calculate counts array 

        Counts are calculated as the gross - background


        Returns
        -------
        counts : np.ndarray
            Array containing the calculated counts

        """
        
        if not len( self.gross ):
            return self.gross.copy()
        else:
            return self.gross - self.background


    @property
    def error(self):
        """ Calculate error array 

        Error is calculated as sqrt( gross + background )

        Returns
        -------
        error : np.ndarray
            Array containing the calculated error

        """

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

        Returns
        -------
        flux_error : np.ndarray
            Array containing the calculated error in flux

        """

        if not len(self.gross):
            return self.gross.copy()
        else:
            return (self.error / self.counts ) *  self.flux


    @property
    def net(self):
        """ Calculate net array 

        Net is calculated as counts / time

        Returns
        -------
        net : np.ndarray
            Array containing calculated net
        
        """



        if not len(self.counts):
            return self.counts.copy()
        else:
            return self.counts / self.times.astype( np.float64 )

        
    @property
    def signal_to_noise(self):
        """ Quick signal to noise estimate

        S/N is calculated as the ratio of gross to error

        Returns
        -------
        sn_ratio : np.ndarray
            Array containing calculated signal to noise ratio

        """

        if not len(self.gross):
            return self.gross.copy()
        else:
            return self.gross / self.error


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

        hdu_out = pyfits.HDUList(pyfits.PrimaryHDU())

        try: hdu_out[0].header = self.hdu[0].header
        except AttributeError: pass 

        hdu_out[0].header['GEN_DATE'] = (str(datetime.now()), 'Creation Date')
        hdu_out[0].header['LC_VER'] = (__version__, 'lightcurve version used')
        hdu_out[0].header['AP_VER'] = (astropy.__version__, 'Astropy version used')
	hdu_out[0].header['NP_VER'] = (np.__version__, 'Numpy version used')
        hdu_out[0].header['SP_VER'] = (scipy.__version__, 'Scipy version used')

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
            
        if outname.endswith('.gz'):
            print("Nope, can't write to gzipped files")
            self.outname = self.outname[:-3]

        hdu_out.writeto( self.outname, clobber=clobber)  

