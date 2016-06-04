"""
Tests for extraction of COS data.

"""
import numpy as np
from astropy.io import fits
import os

from .. import io
from ..cos import extract_index, get_both_filenames

#-------------------------------------------------------------------------------

def generate_test_files(outname='test_corrtag_a.fits', epsilon=1):
    hdu_out = fits.HDUList(fits.PrimaryHDU())
    hdu_out[0].header['detector'] = 'FUV'
    hdu_out[0].header['segment'] = 'FUVA'
    hdu_out[0].header['OBSTYPE'] = 'SPECTROSCOPIC'

    image = np.zeros((1024, 16384))
    image[100] = 1

    y_coords, x_coords = np.where(image > 0)
    n_events = len(y_coords)
    exptime = 100.0

    time_col = fits.Column('time', 'D', 'time', array=np.linspace(0, exptime, n_events))
    rawx_col = fits.Column('rawx', 'I', 'MJD', array=x_coords)
    rawy_col = fits.Column('rawy', 'I', 'MJD', array=y_coords)
    xcorr_col = fits.Column('xcorr', 'I', 'MJD', array=x_coords)
    ycorr_col = fits.Column('ycorr', 'I', 'MJD', array=y_coords)
    xdopp_col = fits.Column('xdopp', 'I', 'counts', array=x_coords)
    xfull_col = fits.Column('xfull', 'I', 'counts', array=x_coords)
    yfull_col = fits.Column('yfull', 'I', 'counts', array=y_coords)
    wavelength_col = fits.Column('wavelength', 'I', 'counts/s', array=np.ones(n_events) * 1200)
    epsilon_col = fits.Column('epsilon', 'D', 'ergs/s', array=np.ones(n_events) * epsilon)
    dq_col = fits.Column('dq', 'I', 'cnts', array=np.zeros(n_events))
    pha_col = fits.Column('pha', 'I', 'counts', array=np.ones(n_events) * 14)

    tab = fits.new_table( [time_col,
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
    hdu_out.append(tab)

    hdu_out[1].header['EXTNAME'] = 'EVENTS'
    hdu_out[1].header['exptime'] = exptime
    hdu_out[1].header['expstart'] = 56000
    for segment in ['a', 'b']:
        hdu_out[1].header['sp_hgt_{}'.format(segment)] = 20
        hdu_out[1].header['sp_loc_{}'.format(segment)] = 100
        hdu_out[1].header['b_hgt1_{}'.format(segment)] = 20
        hdu_out[1].header['b_bkg1_{}'.format(segment)] = 150
        hdu_out[1].header['b_hgt2_{}'.format(segment)] = 20
        hdu_out[1].header['b_bkg2_{}'.format(segment)] = 50

    hdu_out[1].header['sdqflags'] = 16

    hdu_out.writeto(outname, clobber=True)

#-------------------------------------------------------------------------------

class test_extraction:
    def setUp(self):
        self.test_file = 'test_corrtag_a.fits'
        generate_test_files(outname=self.test_file)

    def tearDown(self):
        os.remove(self.test_file)

    def test_int_step(self):
        for stepsize in [.1, .5, 1, 2, 5, 10]:
            lc = io.read(self.test_file, step=stepsize)
            assert lc['gross'].sum() == 16384, "Extraction didn't find all the counts, step={}".format(stepsize)

    def test_step_truncation(self):
        for stepsize in [.3, 3]:
            lc = io.read(self.test_file, step=stepsize)
            assert (lc['gross'].sum() < 16384), "Extraction didn't truncate properly, step={}".format(stepsize)

#-------------------------------------------------------------------------------

class test_epsilon:
    def setUp(self):
        self.test_file = 'test_corrtag_a.fits'
        generate_test_files(outname=self.test_file, epsilon=1.25)

    def tearDown(self):
        os.remove(self.test_file)

    def test_simple_epsilon(self):
        lc = io.read(self.test_file)
        assert lc['gross'].sum() == 16384 * 1.25, 'Espilon not accounted for {}'.format(lc['gross'].sum())

#-------------------------------------------------------------------------------

def test_extract_index():
    """test that the extracted indexes works"""

    test_file = 'dq_corrtag_a.fits'
    generate_test_files(outname=test_file)

    with fits.open(test_file, mode='update') as hdu:
        hdu[1].data['dq'][:] = 8

        good_index = extract_index(hdu, 0, 16384, 0, 1024, -1, 10000, 8)
        assert len( good_index ) == 0, 'Should be no good data: {}'.format( len(good_index) )

        good_index = extract_index( hdu, 0, 16384, 0, 1024, -1, 10000, 8376)
        assert len( good_index ) == 0, 'Still should be no good data: {}'.format( len(good_index) )

    os.remove(test_file)

#-------------------------------------------------------------------------------

def test_get_both_filenames():
    """make sure i can find both"""

    test_names = ('rootname_corrtag_a.fits', 'rootname_corrtag_b.fits')

    assert get_both_filenames(test_names[0]) == test_names, \
        "Failed trying to find the B data with no path"

    assert get_both_filenames(test_names[1]) == test_names, \
        "Failed trying to find the A data with no path"

    test_names = ('/home/ely/rootname_corrtag_a.fits',
                  '/home/ely/rootname_corrtag_b.fits')

    assert get_both_filenames(test_names[0]) == test_names, \
        "Failed trying to find the B data with a path"

    assert get_both_filenames(test_names[1]) == test_names, \
        "Failed trying to find the A data with a path"
