"""General purpose utility functions"""

import os

__all__ = ['expand_refname',
           'enlarge',
           'is_uniq']

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

def enlarge(a, x=2, y=None):
    """Enlarges 2D image array a using simple pixel repetition in both dimensions.
    Enlarges by factor x horizontally and factor y vertically.
    If y is left as None, uses factor x for both dimensions."""

    assert a.ndim == 2

    if y == None:
        y = x

    for factor in (x, y):
        assert factor.__class__ == int
        assert factor > 0

    return a.repeat(y, axis=0).repeat(x, axis=1)

#-------------------------------------------------------------------------------

def is_uniq(values):
    """ Check if input items are unique

    Parameters
    ----------
    values : set
        set of all values

    Returns
    -------
    True/False, MULTI/unique value
    """

    if len(values) == 0:
        return True, ''
    elif len(values) == 1:
        return True, list(values)[0]
    else:
        return False, 'MULTI'

#-------------------------------------------------------------------------------
