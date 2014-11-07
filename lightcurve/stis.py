"""
Holder of Space Telescope Imaging Spectrograph (STIS) classes and utilities

"""

from __future__ import print_function

from astropy.io import fits as pyfits
import os
import numpy as np
from scipy.interpolate import interp1d

from .lightcurve import LightCurve
from .utils import expand_refname

#-------------------------------------------------------------------------------

class StisCurve(LightCurve):
    """
    Subclass for STIS specific routines and abilities

    """

    def __init__(self, **kwargs):
        """ Initialize and extract lightcurve from input corrtag
        """

        if 'step' in kwargs: step = kwargs['step']
        else: step = 1

        if 'filename' in kwargs:
            self.extract(kwargs['filename'], step=step)


    def __str__(self):
        """Prettier representation of object instanct"""

        return "Lightcurve Object of %s" % ( ','.join( self.input_list ) )


    def extract(self, filename, step=1):
        """ Loop over HDUs and extract the lightcurve

        This is the main driver of the lightcuve extracion, and definitely
        needs some better documentation.

        """

        SECOND_PER_MJD = 1.15741e-5

        hdu = pyfits.open(filename)

        time = hdu[1].data['time']
        #time = np.array([round(val, 3) for val in self.hdu[1].data['time']]).astype(np.float64)

        if not len(time):
            end = 0
        else:
            end = min(time.max(),
                      hdu[1].header['EXPTIME'] )


        all_steps = np.arange(0, end+step, step)

        if all_steps[-1] > end:
            truncate = True
        else:
            truncate = False

        gross = np.histogram(time, all_steps)[0]

        self.gross = gross
        self.flux = np.zeros(gross.shape)
        self.background = np.zeros(gross.shape)
        self.mjd = hdu[1].header['EXPSTART'] + np.array(all_steps[:-1]) * SECOND_PER_MJD
        self.bins = np.ones(len(gross)) * step
        self.times = all_steps[:-1]

        if truncate:
            self.gross = self.gross[:-1]
            self.flux = self.flux[:-1]
            self.background = self.background[:-1]
            self.mjd = self.mjd[:-1]
            self.bins = self.bins[:-1]
            self.times = self.times[:-1]

#-------------------------------------------------------------------------------
