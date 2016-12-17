Simple extraction
=================

Performing a basic lightcurve extraction from either COS or STIS individual
datasets is achieved by a simple call to :func:`~lightcurve.read` as shown below.

  >>> import lightcurve
  >>> table = lightcurve.read('lbova4b2q_corrtag_a.fits')

This returns the extracted data as an Astropy Table, which can be used for
further analysis and plotting.

  >>> type(lc)
  astropy.table.table.Table

  >>> lc
  <Table length=750>
  dataset background      mjd       error    bins  ...      flux      signal_to_noise  counts  gross       net
  float64  float32      float64    float32 float64 ...    float64         float32     float32 float32    float64
  ------- ---------- ------------- ------- ------- ... -------------- --------------- ------- ------- -------------
      1.0    12.7314 55899.6041542 46.6337     1.0 ... 0.154727564001         46.3607 2149.24 2161.97 2149.23583984
      1.0    15.3585 55899.6041658 45.9691     1.0 ... 0.149929073704          45.635 2082.44  2097.8 2082.44238281
      1.0     11.119 55899.6041774  45.537     1.0 ... 0.147685074003         45.2928 2051.38  2062.5  2051.3815918
      1.0    16.2978 55899.6041889 45.5294     1.0 ... 0.146895813966         45.1715 2040.33 2056.63 2040.33349609
  ...


Instrument-specific extraction
------------------------------
More specific tailored extractions for each instrument can be seen on their
respective pages below.  These detail specific options for extracion, as well
as how to setup parameters and reference files.

.. toctree::
   :maxdepth: 1

   cos_extraction.rst
   stis_extraction.rst
