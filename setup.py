from distutils.core import setup
import os
import glob

setup(
    name = 'lightcurve',
    url = 'https://github.com/justincely/lightcurve',
    version = '0.0.1', 
    description = 'Create lightcurves from HST/COS data',
    author = 'Justin Ely',
    author_email = 'ely@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    requires = ['numpy','pyfits'],
    scripts =  ['scripts/lightcurve']
    )

