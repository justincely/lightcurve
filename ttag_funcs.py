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
    import numpy as np
   
    SECOND = 1.15741e-5
    hdu = pyfits.open( filename )

    if not len( hdu['events'].data ):
        return ( [], [] )

    DETECTOR = hdu[0].header[ 'DETECTOR' ]
    SEGMENT = hdu[0].header[ 'SEGMENT' ]
    OPT_ELEM = hdu[0].header[ 'OPT_ELEM' ]
    CENWAVE = hdu[0].header[ 'CENWAVE' ]
    APERTURE = hdu[0].header[ 'APERTURE' ]
    SDQFLAGS = hdu[1].header[ 'SDQFLAGS' ]
    
    if not fluxtab:
        fluxtab = hdu[0].header[ 'FLUXTAB' ]

    EXPSTART = hdu[1].header[ 'EXPSTART' ]
    end = hdu['events'].data[ 'time' ].max()

    print "Extracting light curve over", end, 'seconds'

    counts = []
    errors = []
    times = []

    y_start_locs, y_end_locs = get_extraction_regions( hdu )

    print 'Extracting at: ', y_start_locs
    print 'With heights : ', y_end_locs
    steps = range(0, end, step)[:-1]
    N_steps = len( steps )

    if (not xlim) and (DETECTOR == 'FUV'):
        xlim = (0, 16384)
    elif (not xlim) and (DETECTOR == 'NUV'):
        xlim = (0, 1024)
        
    for i,start in enumerate( steps ):
        progress_bar( i, N_steps )
        sub_count = []
        for ystart, yend in zip( y_start_locs, y_end_locs ):
            net = extract_counts(hdu, start, start+step, xlim[0], xlim[1], 
                                 ystart, yend, SDQFLAGS )

            if fluxcal:

                net /=  (float(step) / (xlim[1] - xlim[0]) / height)


                if '$' in fluxtab:
                    fluxpath, fluxfile = fluxtab.split( '$' )
                else:
                    fluxpath = './'
                    fluxfile = fluxtab

                minwave = hdu[1].data[ data_index ]['wavelength'].min()
                maxwave = hdu[1].data[ data_index ]['wavelength'].max()

                #print 'Flux calibrating with: ', fluxfile
                flux_hdu = pyfits.open( fluxfile )
                flux_index = np.where( (flux_hdu[1].data['OPT_ELEM'] == OPT_ELEM) &
                                       (flux_hdu[1].data['CENWAVE'] == CENWAVE) &
                                       (flux_hdu[1].data['APERTURE'] == APERTURE) )[0]

                mean_response_list = []
                for sens_curve,wave_curve in zip( flux_hdu[1].data[flux_index]['SENSITIVITY'],
                                                  flux_hdu[1].data[flux_index]['WAVELENGTH'] ):
                    wave_index = np.where( (wave_curve >= minwave) &
                                           (wave_curve <= maxwave ) )[0]

                    if not len(wave_index): continue

                    mean_response_list.append( sens_curve[ wave_index ].mean() )

                total_fluxcorr = np.mean( mean_response_list )

                net /= total_fluxcorr
                
            sub_count.append( net )

        sample_counts = np.sum( sub_count )

        counts.append( sample_counts )
        errors.append( np.sqrt( sample_counts ) )
        times.append( EXPSTART + (start+1) * SECOND )

    counts = np.array( counts ).astype( np.float64 )
    errors = np.array( errors ).astype( np.float64 )

   
    if normalize:
        errors = errors / counts
        counts = counts / counts.mean()

    return times, counts, errors


def get_extraction_regions( hdu ):
    """Get y_start,y_end for given extraction

    Any extraction regions with values of -999
    will be removed and ignored

    """
    from astroraf.headers import key_exists 


    all_heights = [ hdu[1].header[ 'SP_HGT_%s'% (segment) ] 
                    for segment in ['A','B','C'] if 
                    key_exists( hdu[1].header, 'SP_HGT_%s'% (segment)  ) ]

    all_locations = [ hdu[1].header[ 'SP_LOC_%s'% (segment) ] 
                      for segment in ['A','B','C'] if 
                      key_exists( hdu[1].header,'SP_LOC_%s'% (segment)  ) ]

    for l,h in zip( all_locations, all_heights ):
        if l == -999:
            all_locations.remove( l )
            all_heights.remove( h )

    y_start = [ location - height/2 for location,height in 
                zip( all_locations, all_heights ) ]

    y_end = [ location + height/2 for location,height in 
                zip( all_locations, all_heights ) ]

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
    return len(data_index)
