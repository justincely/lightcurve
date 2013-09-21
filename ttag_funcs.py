"""
Selection of functions dealing with time-tag data in FITS format

"""

from astroraf.misc import progress_bar
import pyfits
import os
import numpy as np

#-------------------------------------------------------------------------------

def ttag_image(in_data, xtype='XCORR', ytype='YCORR', pha=(2, 30), 
               bins=(1024, 16384), times=None, ranges=((0, 1023), (0, 16384)),
               binning=(1, 1), NUV=False):
    """ 
    Bin an events list (*_rawtag_*.fits, *_corrtag_*.fits) into an image.
    
    Events can be filtered out by time, PHA, and/or X and Y ranges
    
    """
    
    from numpy import histogram2d, where, zeros
    from pyfits import getdata

    if isinstance( in_data, str ):
        hdu = pyfits.open( in_data )
    else:
        hdu = in_data

    events = hdu['events'].data

    xlength = (ranges[1][1]+1)
    ylength = (ranges[0][1]+1)
    xbinning = binning[1]
    ybinning = binning[0]

    if NUV:
        bins = (1024, 1024)
        pha = (-1, 1)
        ranges = ( (0, 1023), (0, 1023) )

    if times:
        index = where( (events['TIME'] >= times[0]) & 
                       (events['TIME'] <= times[1]) )
        events =  events[index]

    index = where( (events['PHA'] >= pha[0]) & 
                   (events['PHA'] <= pha[1]) )

    if len(index[0]):
        image = histogram2d( events[ytype][index], 
                             events[xtype][index], 
                             bins=bins, range=ranges)[0]
    else:
        image = zeros( (bins[0]//binning[0], bins[1]//binning[1]) )

    return image

#-------------------------------------------------------------------------------

class lightcurve(object):
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

    def __init__(self, filename, step=5, xlim=None, ylim=None, extract=True,
                 fluxtab=None, normalize=False, writeto=None, clobber=False):

        SECOND = 1.15741e-5
        hdu = pyfits.open( filename )
        self.hdu = hdu
        
        if writeto:
            if isinstance( writeto, bool ):
                self.outname = self.hdu[0].header['rootname'] + '_curve.fits'
            elif isinstance( writeto, str ):
                self.outname = writeto
            else:
                raise NameError( 'writeto should be True/False our string :' )

            if os.path.exists( self.outname ) and not clobber:
                raise IOError( 'Output already exists %s'% self.outname )

        DETECTOR = hdu[0].header[ 'DETECTOR' ]
        OPT_ELEM = hdu[0].header[ 'OPT_ELEM' ]
        CENWAVE = hdu[0].header[ 'CENWAVE' ]
        APERTURE = hdu[0].header[ 'APERTURE' ]
        SDQFLAGS = hdu[1].header[ 'SDQFLAGS' ]

        if DETECTOR == 'NUV':
            if OPT_ELEM == 'G230L':
                hdu_dict = {'A':hdu, 'B':hdu}
            else:
                hdu_dict = {'A':hdu, 'B':hdu, 'C':hdu}

            print 'Making lightcurve for %s'% filename
            print 'Extracting on stripes: ' + ','.join( np.sort(hdu_dict.keys()) )

        elif DETECTOR == 'FUV':
            file_a, file_b = self._get_both_filenames(filename )

            hdu_dict = dict()
            for name, letter in zip( [file_a, file_b], ['A', 'B'] ):
                if os.path.exists( name ):
                    hdu_dict[ letter ] = pyfits.open( name )

            print 'Making lightcurve from: ' + ', '.join( [file_a, file_b] )
            print 'Extracting on stripes: ' +','.join( np.sort(hdu_dict.keys()) )


        print 'DETECTOR: %s'% DETECTOR
        print 'APERTURE: %s'% APERTURE
        print 'OPT_ELEM: %s'% OPT_ELEM
        print 'CENWAVE : %s'% CENWAVE

        if not fluxtab:
            fluxtab = hdu[0].header[ 'FLUXTAB' ]

        EXPSTART = hdu[1].header[ 'EXPSTART' ]
        end = hdu['events'].data[ 'time' ].max()

        print "Extracting over", end, 'seconds'

        all_counts = []
        all_net = []
        all_flux = []
        all_bkgnd = []
        all_error = []
        all_times = []
        all_mjd = []

        all_times = np.array( range(0, end, step)[:-1] )
        N_steps = len( all_times )

        if (not xlim) and (DETECTOR == 'FUV'):
            xlim = (0, 16384)
        elif (not xlim) and (DETECTOR == 'NUV'):
            xlim = (0, 1024)

        if (not ylim) and (DETECTOR == 'FUV'):
            ylim = 1024
        elif (not ylim) and (DETECTOR == 'NUV'):
            ylim = 512

        for i,start in enumerate( all_times ):
            progress_bar( i, N_steps )
            sub_count = []
            sub_bkgnd = []
            sub_net = []
            sub_flux = []
            for segment, hdu in hdu_dict.iteritems():
                ystart, yend = self._get_extraction_region( hdu, segment, 'spectrum' )
                counts, wmin, wmax = self._extract_counts(hdu, start, 
                                                          start+step, 
                                                          xlim[0], xlim[1], 
                                                          ystart, yend, 
                                                          SDQFLAGS )

                bstart, bend = self._get_extraction_region( hdu, segment, 'background1' )
                b_counts, b_wmin, b_wmax = self._extract_counts(hdu, start, 
                                                                start+step, 
                                                                xlim[0], xlim[1], 
                                                                bstart, bend, 
                                                                SDQFLAGS )

                b_correction = ( (bend - bstart) / (yend -ystart) )

                background = b_counts * b_correction

                net = float(counts) / (yend-ystart) 
                flux = net / self._get_flux_correction( fluxtab, OPT_ELEM, 
                                                        CENWAVE, APERTURE, 
                                                        wmin, wmax )

                sub_count.append( counts )
                sub_bkgnd.append( background )
                sub_net.append( net )
                sub_flux.append( flux )

            sample_counts = np.sum( sub_count )
            sample_bkgnd = np.sum( sub_bkgnd )

            all_counts.append( sample_counts - sample_bkgnd )
            all_net.append( np.sum(sub_net) )
            all_flux.append( np.sum(sub_flux) )
            all_bkgnd.append( sample_bkgnd )
            all_error.append( np.sqrt( sample_counts + sample_bkgnd ) )
            all_mjd.append( EXPSTART + (start * SECOND) )

        self.counts = np.array( all_counts ).astype( np.float64 )
        self.net = np.array( all_net ).astype( np.float64 )
        self.flux = np.array( all_flux ).astype( np.float64 )
        self.background = np.array( all_bkgnd ).astype( np.float64 )
        self.mjd = np.array( all_mjd ).astype( np.float64 )
        self.times = np.array( all_times ).astype( np.float64 )
        self.error = np.array( all_error ).astype( np.float64 )


        if normalize:
            self.error = self.error / self.counts
            self.counts = self.counts / self.counts.mean()
            self.net = self.net / self.net.mean()
            self.background = self.background / self.background.mean()
            self.flux = self.flux / self.flux.mean()

        if writeto:
            self.write( clobber=clobber  )


    def __add__( self, other ):
        self.counts = np.concatenate( [self.counts, other.counts] )
        self.net = np.concatenate( [self.net, other.counts] )
        self.flux = np.concatenate( [self.flux, other.flux] )
        self.background = np.concatenate( [self.background, other.background] )
        self.mjd = np.concatenate( [self.mjd, other.mjd] )
        self.times = np.concatenate( [self.times, other.times + self.times[-1]] )
        self.error = np.concatenate( [self.error, other.error] )

        return self


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
            height =  hdu[1].header['SP_HGT_%s'% (segment)]
            location = hdu[1].header['SP_LOC_%s'% (segment)]

        elif mode == 'background1':
            height =  hdu[1].header['B_HGT1_%s'% (segment)]
            location = hdu[1].header['B_BKG1_%s'% (segment)] 

        elif mode == 'background2':
            height =  hdu[1].header['B_HGT2_%s'% (segment)]
            location = hdu[1].header['B_BKG2_%s'% (segment)] 

        y_start = location - height/2
        y_end = location + height/2

        return y_start, y_end


    def _extract_counts(self, hdu, start, end, x_start, x_end, 
                        y_start, y_end, sdqflags=0):
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


    def _get_mean_pha(self, hdu, start, end, x_start, x_end, 
                      y_start, y_end, sdqflags=0):
        """
        Find mean PHA of events

        """

        data_index = np.where( ( hdu[1].data['TIME'] >= start ) & 
                               ( hdu[1].data['TIME'] < end ) &
                               ( hdu[1].data['XCORR'] >= x_start ) & 
                               ( hdu[1].data['XCORR'] < x_end ) &
                               ( hdu[1].data['YCORR'] >= y_start ) & 
                               ( hdu[1].data['YCORR'] < y_end ) &
                               ~( hdu[1].data['DQ'] & sdqflags ) &
                               ( (hdu[1].data['WAVELENGTH'] > 1217) | 
                                 (hdu[1].data['WAVELENGTH'] < 1214) ) )[0]

        return hdu[1].data['PHA'][index].mean()


    def _get_flux_correction(self, fluxtab, opt_elem, cenwave, 
                             aperture, minwave, maxwave):
        """ Return integrated flux calibration over given wavelength range 
        """

        if '$' in fluxtab:
            fluxpath, fluxfile = fluxtab.split( '$' )
        else:
            fluxpath = './'
            fluxfile = fluxtab

        if not os.path.exists( fluxfile ):
            #print 'WARNING: Fluxfile not available, using unity flux calibration'
            return 1

        flux_hdu = pyfits.open( fluxfile )
        flux_index = np.where( (flux_hdu[1].data['OPT_ELEM'] == opt_elem) &
                               (flux_hdu[1].data['CENWAVE'] == cenwave) &
                               (flux_hdu[1].data['APERTURE'] == aperture) )[0]

        mean_response_list = []
        for sens_curve,wave_curve in zip( flux_hdu[1].data[flux_index]['SENSITIVITY'],
                                          flux_hdu[1].data[flux_index]['WAVELENGTH'] ):
            wave_index = np.where( (wave_curve >= minwave) &
                                   (wave_curve <= maxwave ) )[0]

            if not len(wave_index): continue

            mean_response_list.append( sens_curve[ wave_index ].mean() )

        total_fluxcorr = np.mean( mean_response_list )

        return total_fluxcorr


    def write(self, outname=None, clobber=False):
        """ Write out to FITS file
        """

        if isinstance( outname, str ):
            self.outname = outname


        hdu_out = pyfits.HDUList(pyfits.PrimaryHDU())

        keyword_list = ['TELESCOP', 'INSTRUME', 'DETECTOR', 'OPT_ELEM',
                        'ROOTNAME', 'TARGNAME', 'RA_TARG', 'DEC_TARG',
                        'PROPOSID', 'OPUS_VER', 'CAL_VER', 'PROCTIME',
                        'APERTURE', 'OPT_ELEM', 'CENWAVE', ]

        for kw in keyword_list:
            hdu_out[0].header[kw] = self.hdu[0].header[kw]
 
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

        hdu_out.writeto( self.outname, clobber=clobber )  
