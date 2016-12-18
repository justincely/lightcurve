from scipy.signal import lombscargle

__all__ = ['lomb']

#-------------------------------------------------------------------------------

def lomb(time, counts, frequencies):
    """Compute the lombscargle periodigram

    Necessary wrapper around the set lomscargle algorithm
    https://github.com/scipy/scipy/issues/2643

    Parameters
    ----------
    time : np.ndarray
        array of data times
    counts : np.ndarray
        array of counts
    frequencies : np.ndarray
        What frequencies

    Returns
    -------
    freqs : np.ndarray
        calculated freqencies

    """

    time = time.byteswap().newbyteorder().astype('float64')
    counts = counts.byteswap().newbyteorder().astype('float64')

    freqs = lombscargle(time, counts, frequencies)

    return freqs

#-------------------------------------------------------------------------------
