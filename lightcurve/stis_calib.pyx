import numpy as np
cimport numpy as np

def map_image(np.ndarray image, np.ndarray xcoords, np.ndarray ycoords):
    print("Starting the Cython stuff")
    cdef int x=0, y=0, n_coord=0

    xdim = xcoords.shape[0]
    ydim = xcoords.shape[1]
    n_coord = len(xcoords)
    out_vals = np.zeros(n_coord)

    for i in range(n_coord):
        x = xcoords[i]
        y = ycoords[i]
        if (not 0 < x < 2047) or (not 0 < y < 2047):
            out_vals[i] = 0
        else:
            out_vals[i] = image[x, y]

    return out_vals
