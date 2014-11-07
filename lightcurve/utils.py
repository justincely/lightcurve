#-------------------------------------------------------------------------------

def expand_refname(refname):
    '''Expand header reference file name to full path if $ is present.

    Parameters
    ----------
    refname, str
        reference file name

    Returns
    -------
    reffile, str
        expanded full path to reference file

    '''

    if '$' in refname:
        refpath, reffile = refname.split('$')

        try:
            reffile = os.path.join(os.environ[refpath], reffile)
        except KeyError:
            pass

    else:
        refpath = './'
        reffile = refname

    return reffile

#-------------------------------------------------------------------------------
