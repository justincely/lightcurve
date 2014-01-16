"""
Holder of Cosmic Origins Spectrograph (COS) classes and utilities

"""

from __future__ import print_function

__all__ = ['extract_index']

from astropy.io import fits as pyfits
import os
import numpy as np
from scipy.interpolate import interp1d

from .lightcurve import LightCurve

#-------------------------------------------------------------------------------

class CosCurve( LightCurve ):
    """
    Subclass for COS specific routines and abilities

    """
      
    def __init__(self, **kwargs):
        """ Initialize and extract lightcurve from input corrtag
        """

        if 'step' in kwargs: step = kwargs['step']
        else: step = 1

        if 'wlim' in kwargs: wlim = kwargs['wlim']
        else: wlim = (1, 10000)

        if 'filename' in kwargs:
            self.read_cos( kwargs['filename'] )
            self.extract( step=step, wlim=wlim )
        else:
            pass


    def __str__(self):
        """Prettier representation of object instanct"""

        return "Lightcurve Object of %s" % ( ','.join( self.input_list ) ) 


    def read_cos(self, filename ):
        """ Extract light curve from COS data"""

        self.hdu = pyfits.open( filename )
        self.input_filename = filename

        self._get_hdus()


    def extract(self, step=1, xlim=(0, 16384), wlim=(2, 10000), ylim=None ):
        """ Loop over HDUs and extract the lightcurve

        """

        SECOND_PER_MJD = 1.15741e-5


        if not len( self.hdu[1].data['time'] ): 
            end = 0
        else:
            end = min( self.hdu['events'].data[ 'time' ].max(), 
                       self.hdu[1].header['EXPTIME'] )


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
            #print('Source fluxcal')
            if not ylim:
                ystart, yend = self._get_extraction_region( hdu, segment, 'spectrum' )
            else:
                ystart, yend = ylim[0], ylim[1]
            #print(segment, ystart, yend)
            #print(xlim, ystart, yend, wlim, hdu[1].header['sdqflags'])
            index = extract_index(hdu,
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
            #print('Background fluxcal')
            bstart, bend = self._get_extraction_region( hdu, segment, 'background1' )
            index = extract_index(hdu,
                                  xlim[0], xlim[1], 
                                  bstart, bend, 
                                  wlim[0], wlim[1],
                                  hdu[1].header['sdqflags'] )

            bstart, bend = self._get_extraction_region( hdu, segment, 'background2' )
            index = np.hstack( (index, extract_index(hdu,
                                                     xlim[0], xlim[1], 
                                                     bstart, bend, 
                                                     wlim[0], wlim[1],
                                                     hdu[1].header['sdqflags']) ) )
            
            b_corr = ( (bend - bstart) / (yend -ystart) ) / 2.
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


    def _get_tds(self, hdu, index):
        """ return tds correction for each event
        """

        try:
            tdstab = hdu[0].header['TDSTAB']
        except KeyError:
            tdstab = ''

            
        tdsfile = expand_refname( tdstab )

        if (not tdsfile) or (not os.path.exists( tdsfile ) ):
            print(' WARNING: tdsfile not available %s,' % tdsfile )
            print(' using unity TDS correction instead.' )
            return np.ones( hdu[ 'events' ].data['time'].shape )[index]
        
        tds_data = pyfits.getdata( tdsfile, ext=1 )
        REF_TIME = pyfits.getval( tdsfile, 'REF_TIME', ext=1)
        mode_index = np.where( (tds_data['OPT_ELEM'] == hdu[0].header['opt_elem'] ) &
                               (tds_data['APERTURE'] == hdu[0].header['aperture'] ) &
                               (tds_data['SEGMENT'] == hdu[0].header['segment'] ) )[0]


        if len( mode_index ) == 0:
            print('No row in tdstab found for this dataset, no TDSCORR performed')
            return np.ones( hdu[ 'events' ].data['time'].shape )[index]
        elif len( mode_index ) > 1:
            raise ValueError('Too many rows found: {}'.format( len(mode_index) ) )

        mode_line = tds_data[ mode_index ]

        tds_nt = mode_line['NT'][0]
        tds_wavelength = mode_line['WAVELENGTH'][0]
        smaller_index = np.where( mode_line['TIME'][0][:tds_nt] < hdu[1].header['EXPSTART'] )[0]
        time_index = smaller_index.max()

        tds_slope = mode_line['SLOPE'][0]
        tds_intercept = mode_line['INTERCEPT'][0]

        calculate_drop = lambda time, ref_time, slope, intercept: \
            ( (( time - ref_time ) * slope) / (365.25*100) ) + intercept

        correction_array = np.zeros( len(tds_wavelength) )
        for i, wave in enumerate( tds_wavelength ):
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


        fluxfile = expand_refname( fluxtab )

        if (not fluxfile) or (not os.path.exists( fluxfile ) ):
            print(' WARNING: Fluxfile not available %s,' % fluxfile )
            print(' using unity flux calibration instead.' )
            return np.ones( hdu[ 'events' ].data['time'].shape )[index]


        flux_hdu = pyfits.open( fluxfile )
        setting_index = np.where( (flux_hdu[1].data['SEGMENT'] == hdu[0].header['segment']) &
                               (flux_hdu[1].data['OPT_ELEM'] == hdu[0].header['opt_elem']) &
                               (flux_hdu[1].data['CENWAVE'] == hdu[0].header['cenwave']) &
                               (flux_hdu[1].data['APERTURE'] == hdu[0].header['aperture']) )[0]
        if len( setting_index ) == 0:
            print('No row in fluxtab found for this dataset, no FLUXCAL performed')
            return np.ones( hdu[ 'events' ].data['time'].shape )[index]
        elif len( setting_index ) > 1:
            raise ValueError('Too many rows found: {}'.format( len(setting_index) ) )
        

        resp_wave = flux_hdu[1].data[setting_index]['WAVELENGTH'].flatten()
        response = flux_hdu[1].data[setting_index]['SENSITIVITY'].flatten()

        data_max = hdu['events'].data['wavelength'][index].max()
        data_min = hdu['events'].data['wavelength'][index].min()

        if data_max > resp_wave.max():
            print( "Expanding minumum response curve by {}".format( data_max - resp_wave.max() ) )
            resp_wave[ np.argmax( resp_wave ) ] = data_max

        if data_min < resp_wave.min():
            print( "Expanding minimum response curve by {}".format( data_min - resp_wave.min() ) )
            resp_wave[ np.argmin( resp_wave ) ] = data_min

        interp_func = interp1d( resp_wave, response )
        all_resp = interp_func( hdu[ 'events' ].data['wavelength'][index] )

        return all_resp

#--------------------------------------------------------------

def extract_index( hdu, x_start, x_end, 
                   y_start, y_end, w_start, w_end, sdqflags=0):
    """
    Extract event indeces from given HDU using input parameters.

    Wavelength regions containing geocoronal airglow emission will be excluded
    automatically.  These regions include Lyman Alpha between 1214 and 1217
    Angstroms and Oxygen I from 1300 to 1308.


    Parameters
    ----------
    hdu : HDUlist
        Header Data Unit from COS corrtag file
        
    x_start : float
        Lower bound of events in pixel space

    x_end : float
        Upper bound of events in pixel space

    y_start : float
        Lower bound of events in pixel space

    y_end : float
        Upper bound of events in pixel space

    w_start : float
        Lower bound of events in wavelength space

    w_end : float
        Upper bound of events in wavelength space

    sdqflags : int
        Bitwise DQ value of bad events

    Returns
    -------
    data_index : np.ndarray
        Indeces of desired events

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
                             (hdu[1].data['WAVELENGTH'] < 1300) ) 
                           ) [0]

    return data_index

#--------------------------------------------------------------

def expand_refname( refname ):
    '''  Expand header reference file name to full path if $ is
    present.

    '''

    if '$' in refname:
        refpath, reffile = refname.split( '$' )
        
        try:
            reffile = os.path.join( os.environ[refpath], reffile )
        except KeyError:
            pass

    else:
        refpath = './'
        reffile = refname

    return reffile

