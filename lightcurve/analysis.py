from scipy.signal import lombscargle

#-------------------------------------------------------------------------------

def lomb(time, counts, frequencies):
    """Compute the lombscargle periodigram

    Necessary wrapper around the set lomscargle algorithm
    https://github.com/scipy/scipy/issues/2643
    """

    time = time.byteswap().newbyteorder().astype('float64')
    counts = counts.byteswap().newbyteorder().astype('float64')

    freqs = lombscargle(time, counts, frequencies)

    return freqs

#-------------------------------------------------------------------------------
