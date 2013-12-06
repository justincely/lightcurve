"""
simplify IO

"""

from astropy.io import fits as pyfits

from .lightcurve import LightCurve
from .cos import CosCurve

__all__ = ['open']

#--------------------------------------------------------------

def open( **kwargs ):
    """ Open file into lightcurve
    
    """

    if not 'filename' in kwargs:
        raise IOError( 'filename must be supplied' )

    filetype = check_filetype( kwargs['filename'] )

    if filetype == 'corrtag':
        return CosCurve( **kwargs )
    elif filetype == 'lightcurve':
        return open_lightcurve( kwargs['filename'] )

#--------------------------------------------------------------

def check_filetype(filename):
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

#--------------------------------------------------------------

def open_lightcurve(filename):    
    """ Read lightcurve from fits file back into base object"""

    out_obj = LightCurve()

    hdu = pyfits.open( filename)

    out_obj.times = hdu[1].data['times']
    out_obj.gross = hdu[1].data['gross']
    out_obj.mjd = hdu[1].data['mjd']
    out_obj.flux = hdu[1].data['flux']
    out_obj.background = hdu[1].data['background']

    return out_obj
