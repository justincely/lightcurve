"""
Tests for extraction of COS data.  

"""

from ..lightcurve import LightCurve
import numpy as np
import pyfits
import os

#-------------------------------------------------------------------------------

def generate_test_files():
    outname = 'test_corrtag_a.fits'

    hdu_out = pyfits.HDUList(pyfits.PrimaryHDU())
    hdu_out[0].header.update( 'detector', 'FUV' )
    hdu_out[0].header.update( 'segment', 'FUVA' )

    image = np.zeros( (1024, 16384) )
    image[100] = 1

    y_coords, x_coords = np.where( image > 0 )
    n_events = len( y_coords )
    exptime = 100.0

    time_col = pyfits.Column('time', 'D', 'time', array=np.linspace(0, exptime, n_events) )
    rawx_col = pyfits.Column('rawx', 'I', 'MJD', array=x_coords ) 
    rawy_col = pyfits.Column('rawy', 'I', 'MJD', array=y_coords ) 
    xcorr_col = pyfits.Column('xcorr', 'I', 'MJD', array=x_coords ) 
    ycorr_col = pyfits.Column('ycorr', 'I', 'MJD', array=y_coords ) 
    xdopp_col = pyfits.Column('xdopp', 'I', 'counts', array=x_coords )    
    xfull_col = pyfits.Column('xfull', 'I', 'counts', array=x_coords )
    yfull_col = pyfits.Column('yfull', 'I', 'counts', array=y_coords )
    wavelength_col = pyfits.Column('wavelength', 'I', 'counts/s', array=np.ones( n_events ) * 1400 )
    epsilon_col = pyfits.Column('epsilon', 'D', 'ergs/s', array=np.ones( n_events ) )
    dq_col = pyfits.Column('dq', 'I', 'cnts', array=np.zeros( n_events ) )
    pha_col = pyfits.Column('pha', 'I', 'counts', array=np.ones( n_events ) * 14 )
    
    tab = pyfits.new_table( [time_col,
                             rawx_col,
                             rawy_col,
                             xcorr_col,
                             ycorr_col,
                             xdopp_col,
                             xfull_col,
                             yfull_col,
                             wavelength_col,
                             epsilon_col,
                             dq_col,
                             pha_col] )
    hdu_out.append( tab )

    hdu_out[1].header.update( 'EXTNAME', 'EVENTS' )
    hdu_out[1].header.update( 'exptime', exptime )
    hdu_out[1].header.update( 'expstart', 56000 )
    for segment in ['a', 'b']:
        hdu_out[1].header.update( 'sp_hgt_{}'.format( segment ), 20 )
        hdu_out[1].header.update( 'sp_loc_{}'.format( segment ), 100 )
        hdu_out[1].header.update( 'b_hgt1_{}'.format( segment ), 20 )
        hdu_out[1].header.update( 'b_bkg1_{}'.format( segment ), 150 )

    hdu_out[1].header.update( 'sdqflags', 16 )

    hdu_out.writeto( outname, clobber=True )  

#-------------------------------------------------------------------------------

def test_FUV():
    """ Test the FUV file extraction """

    generate_test_files()
    lc = LightCurve(filename='test_corrtag_a.fits')

    for stepsize in [.1, .5, 1, 2, 5, 10]:
        lc.extract( step=stepsize )
        assert lc.gross.sum() == 16384, "Extraction didn't find all the counts, step={}".format(step)

    ### uneven steps, make sure last bit is truncated.
    ### This should also check that the are above some value too
    for stepsize in [.3, 3]:
        lc.extract( step=stepsize )
        assert (lc.gross.sum() < 16384), "Extraction didn't truncate properly, step={}".format(step)
