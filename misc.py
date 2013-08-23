__all__ = [ 'remake_asn', 'read_asn' ]

import pyfits
import calcos

def remake_asn( asn_name, member_ext='_x1d.fits', product_ext='_x1dsum.fits',
                allow_missing=False ):
    """
    Use the fpavg function to re-create the x1dsum

    """

    from calcos import fpavg
    import os

    if not 'lref' in os.environ.keys():
        os.environ['lref'] = '/grp/hst/cdbs/lref/'

    members, product = read_asn( asn_name )

    missing_members = [ item for item in members if os.path.exists( item ) ]

    if len( missing_members ) and (not allow_missing):
        raise IOError( "The following members were not found\n %s"% 
                       (','.join( missing_members) ) )

    elif len( missing_members ) and allow_missing:
        members = [ item for item in members if not (item in missing_members) ]


    if len( product ) > 2:
        raise IOError( 'Too many products' )
    else:
        product = product[0] + product_ext

    members = [ item + member_ext for item in members ]

    fpavg.fpAvgSpec( members, product)

    
def read_asn( asn_name ):
    """

    Reads in an input association table and returns a tuple
    of the members and products.
    
    Inputs:
        asn
        either string or open fits file

    Output:
        (list of members,list of products)
    """


    asn_data = pyfits.open( asn_name )

    members = [ line['MEMNAME'].lower() for line in asn_data['ASN'].data 
                if line['MEMTYPE'] == 'EXP-FP' ]
    products = [ line['MEMNAME'].lower() for line in asn_data['ASN'].data 
                 if line['MEMTYPE'] == 'PROD-FP' ]

    return members, products
