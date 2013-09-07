"""
Selection of functions dealing with time-tag data in FITS format

"""

#-------------------------------------------------------------------------------

def ttag_image(in_data, xtype='XCORR', ytype='YCORR', pha=(2, 30), 
               bins=(1024, 16384), times=None, ranges=((0, 1023), (0, 16384)),
               binning=(1, 1), NUV=False):
    """ 
    Bin an events list (*_rawtag_*.fits, *_corrtag_*.fits) into an image.
    
    Events can be filtered out by time, PHA, and/or X and Y ranges
    
    """
    
    try: histogram2d
    except NameError: from numpy import histogram2d, where, zeros
    
    try: getdata
    except NameError: from pyfits import getdata
    
    try: events = getdata(in_data, 1)
    except: events = in_data

    xlength = (ranges[1][1]+1)
    ylength = (ranges[0][1]+1)
    xbinning = binning[1]
    ybinning = binning[0]

    if NUV:
        bins = (1024, 1024)
        pha = (-1, 1)
        ranges = ( (0, 1023), (0, 1023) )

    if times:
        index = where( (events['TIME']>=times[0]) & 
                       (events['TIME'] <= times[1]) )
        events =  events[index]

    index = where((events['PHA']>=pha[0])&(events['PHA']<=pha[1]))

    if len(index[0]):
        image = histogram2d( events[ytype][index], 
                                       events[xtype][index], 
                                       bins=bins, range=ranges)[0]
    else:
        image = zeros( (bins[0]//binning[0], bins[1]//binning[1]) )

    return image

#-------------------------------------------------------------------------------

def lightcurve( filename, step=5, xlim=None, ylim=None, extract=True,
                fluxcal=False, fluxtab=None, normalize=False):
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
    from astroraf.misc import progress_bar
    import pyfits
    import os
    import numpy as np
   
    SECOND = 1.15741e-5
    hdu = pyfits.open( filename )

    DETECTOR = hdu[0].header[ 'DETECTOR' ]

    if DETECTOR == 'NUV':
        hdu_dict = {'A':hdu, 'B':hdu, 'C':hdu}
    elif DETECTOR == 'FUV':
        file_a, file_b = get_both_filenames( filename )

        hdu_dict = dict()
        for name, letter in zip( [file_a, file_b], ['A', 'B'] ):
            if os.path.exists( name ):
                print name
                hdu_dict[ letter ] = pyfits.open( name )

    OPT_ELEM = hdu[0].header[ 'OPT_ELEM' ]
    CENWAVE = hdu[0].header[ 'CENWAVE' ]
    APERTURE = hdu[0].header[ 'APERTURE' ]
    SDQFLAGS = hdu[1].header[ 'SDQFLAGS' ]
    
    if not fluxtab:
        fluxtab = hdu[0].header[ 'FLUXTAB' ]

    EXPSTART = hdu[1].header[ 'EXPSTART' ]
    end = hdu['events'].data[ 'time' ].max()

    print "Extracting light curve over", end, 'seconds'

    all_counts = []
    all_net = []
    all_flux = []
    all_errors = []
    all_times = []

    steps = range(0, end, step)[:-1]
    N_steps = len( steps )

    if (not xlim) and (DETECTOR == 'FUV'):
        xlim = (0, 16384)
    elif (not xlim) and (DETECTOR == 'NUV'):
        xlim = (0, 1024)

    if (not ylim) and (DETECTOR == 'FUV'):
        ylim = 1024
    elif (not ylim) and (DETECTOR == 'NUV'):
        ylim = 512
        
    for i,start in enumerate( steps ):
        progress_bar( i, N_steps )
        sub_count = []
        sub_net = []
        sub_flux = []
        for segment, hdu in hdu_dict.iteritems():
            ystart, yend = get_extraction_region( hdu, segment )
            counts, wmin, wmax = extract_counts(hdu, start, start+step, xlim[0], xlim[1], 
                                             ystart, yend, SDQFLAGS )

            net = float(counts) / (yend-ystart) 
            flux = net / get_flux_correction( fluxtab, OPT_ELEM, 
                                              CENWAVE, APERTURE, wmin, wmax )

            sub_count.append( counts )
            sub_net.append( net )
            sub_flux.append( flux )

        sample_counts = np.sum( sub_count )

        all_counts.append( sample_counts )
        all_net.append( np.sum(sub_net) )
        all_flux.append( np.sum(sub_flux) )
        all_errors.append( np.sqrt( sample_counts ) )
        all_times.append( EXPSTART + (start+1) * SECOND )

    counts = np.array( all_counts ).astype( np.float64 )
    net = np.array( all_net ).astype( np.float64 )
    flux = np.array( all_flux ).astype( np.float64 )
    times = np.array( all_times ).astype( np.float64 )
    errors = np.array( all_errors ).astype( np.float64 )

   
    if normalize:
        errors = errors / counts
        counts = counts / counts.mean()
        net = net / net.mean()
        flux = flux / flux.mean()

    return times, counts, errors


def get_both_filenames( filename ):
    import os
    import pyfits

    assert pyfits.getval( filename, 'DETECTOR' ) == 'FUV', 'This only works for FUV data'

    if pyfits.getval(filename, 'SEGMENT') == 'FUVA':
        other_filename = filename.replace('_a.fits', '_b.fits')
    elif pyfits.getval(filename, 'SEGMENT') == 'FUVB':
        other_filename = filename.replace('_b.fits', '_a.fits')
        
    filename_list = [filename, other_filename]
    filename_list.sort()

    return ( filename_list[0], filename_list[1] )


def get_extraction_region( hdu, segment ):
    """Get y_start,y_end for given extraction

    """

    height =  hdu[1].header['SP_HGT_%s'% (segment)]
    location = hdu[1].header['SP_LOC_%s'% (segment)]

    y_start = location - height/2
    y_end = location + height/2

    return y_start, y_end


def extract_counts(hdu, start, end, x_start, x_end, y_start, y_end, sdqflags=0):
    """
    Extract counts from given HDU using input parameters.

    Lyman Alpha line counts will be excluded by default.

    """
    import numpy as np
    data_index = np.where( ( hdu[1].data['TIME'] >= start ) & 
                           ( hdu[1].data['TIME'] < end ) &
                           ( hdu[1].data['XCORR'] > x_start ) & 
                           ( hdu[1].data['XCORR'] < x_end ) &
                           ( hdu[1].data['YCORR'] > y_start ) & 
                           ( hdu[1].data['YCORR'] < y_end ) &
                           ~( hdu[1].data['DQ'] & sdqflags ) &
                           ( (hdu[1].data['WAVELENGTH'] > 1217) | 
                             (hdu[1].data['WAVELENGTH'] < 1214) ) )[0]

    minwave = hdu[1].data[ data_index ]['wavelength'].min()
    maxwave = hdu[1].data[ data_index ]['wavelength'].max()

    return len(data_index), minwave, maxwave


def get_flux_correction( fluxtab, opt_elem, cenwave, aperture, minwave, maxwave):
    import pyfits
    import numpy as np
    import os

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
