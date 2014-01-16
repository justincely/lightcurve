"""
I/O operations to get data into a LightCurve

"""

from astropy.io import fits as pyfits

from .lightcurve import LightCurve
from .cos import CosCurve

__all__ = ['open']

#--------------------------------------------------------------

def open( **kwargs ):
    """ Open file into lightcurve

    filename must be supplied in kwargs

    Parameters
    ----------
    **kwargs : 
        Additional arguements to be passed to lightcurve instantiations

    Returns
    -------
    LightCurve or subclass

    Raises
    ------
    IOError
        If filename is not supplied in **kwargs

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

    lightcurve_names = set( ['TIMES',
                             'MJD',
                             'GROSS',
                             'COUNTS',
                             'NET',
                             'FLUX',
                             'FLUX_ERROR',
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
    """ Read lightcurve from fits file back into base object

    Parameters
    ----------
    filename : str
        Filename of FITS lightcurve to open

    Returns
    -------
    out_lc : LightCurve
        LightCurve instantiation containing data from file

    """

    out_lc = LightCurve()

    hdu = pyfits.open( filename)

    out_lc.times = hdu[1].data['times']
    out_lc.gross = hdu[1].data['gross']
    out_lc.mjd = hdu[1].data['mjd']
    out_lc.flux = hdu[1].data['flux']
    out_lc.background = hdu[1].data['background']

    return out_obj
