from distutils.core import setup, Extension
import os
import glob
import numpy as np

stis_cal_module = Extension('lightcurve/stis_calib',
                            sources=['lightcurve/stis_calib.pyx'])

setup(
    name = 'lightcurve',
    url = 'http://justincely.github.io/lightcurve/',
    version = '0.5.1',
    description = 'Create lightcurves from HST/COS data',
    author = 'Justin Ely',
    author_email = 'ely@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    packages = ['lightcurve'],
    requires = ['numpy','scipy','astropy', 'cython'],
    scripts =  ['scripts/lightcurve'],
    include_dirs =[np.get_include()],
    )
