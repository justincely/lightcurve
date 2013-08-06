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

def lightcurve( filename, xlim=(0, 1024), ylim=(0, 1024), step=5, extract=True,
                fluxcal=False, fluxtab=None, normalize=False):
    """
    Turn an event list (*_rawtag_*.fits, *_corrtag_*.fits) into a lightcurve
    
    """
    from astroraf.headers import key_exists 
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
    
    if not fluxtab:
        fluxtab = hdu[0].header[ 'FLUXTAB' ]

    start = hdu[1].header[ 'EXPSTART' ]
    end = hdu['events'].data[ 'time' ].max()

    print "Extracting light curve over", end, 'seconds'

    counts = []
    errors = []
    times = []

    all_heights = [ hdu[1].header[ 'SP_HGT_%s'% (segment) ] 
                    for segment in ['A','B','C'] if 
                    key_exists( hdu[1].header, 'SP_HGT_%s'% (segment)  ) ]
  
    all_locations = [ hdu[1].header[ 'SP_LOC_%s'% (segment) ] 
                      for segment in ['A','B','C'] if 
                      key_exists( hdu[1].header,'SP_LOC_%s'% (segment)  ) ]

    print 'Extracting at: ', all_locations
    for i in range(0, end, step)[:-1]:
        print i
        sub_count = 0
        for loc, height in zip( all_locations, all_heights ):
            index = np.where( (hdu[1].data['TIME'] >= i) & 
                              ( hdu[1].data['TIME'] < i+step) &
                              (hdu[1].data['XCORR'] > xlim[0]) & 
                              ( hdu[1].data['XCORR'] < xlim[1]) &
                              (hdu[1].data['YCORR'] > loc-(height/2) ) & 
                              ( hdu[1].data['YCORR'] < loc+(height/2) ) &
                              ( (hdu[1].data['WAVELENGTH'] > 1217) | 
                                (hdu[1].data['WAVELENGTH'] < 1214) ) )[0]

            sub_count += len(index) / float(step) / (xlim[1] - xlim[0]) / height
            print 100 * np.sqrt(len(index) ) / len(index)

        counts.append( sub_count )
        errors.append( 100 * np.sqrt( len(index) ) / len(index) )
        times.append( start + (i+1) * SECOND )

    counts = np.array( counts )
    errors = np.array( errors )

    if fluxcal:
        if '$' in fluxtab:
            fluxpath, fluxfile = fluxtab.split( '$' )
        else:
            fluxpath = './'
            fluxfile = fluxtab

        print 'Flux calibrating with: ', fluxfile
        flux_hdu = pyfits.open( fluxfile )
        index = np.where( (flux_hdu[1].data['OPT_ELEM'] == OPT_ELEM) &
                          (flux_hdu[1].data['CENWAVE'] == CENWAVE) &
                          (flux_hdu[1].data['APERTURE'] == APERTURE) )[0]

        total_fluxcorr = np.mean( [ flux_hdu[1].data[i]['SENSITIVITY'] 
                                    for i in index ] )
        counts /= total_fluxcorr            
        errors /= total_fluxcorr

    if normalize:
        counts /= np.median( counts )

    return times, counts, errors
