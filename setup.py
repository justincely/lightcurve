from distutils.core import setup, Extension
import os
import glob

stis_cal_module = Extension('lightcurve/stis_cal',
                            sources = ['lightcurve/stis_cal.cc'],
                            extra_link_args=[ '-pthread',
                                             '-shared',
                                             '-lm']
                             )

setup(
    name = 'lightcurve',
    url = 'http://justincely.github.io/lightcurve/',
    version = '0.2.1',
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
    packages = ['lightcurve'],
    requires = ['numpy','scipy','astropy'],
    scripts =  ['scripts/lightcurve'],
    ext_modules = [stis_cal_module]
    )
