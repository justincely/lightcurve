cimport numpy as np

def map_image(image, xcoords, ycoords):
    print("Starting the Cython stuff")
    out_vals = np.zeros(xcoords.shape)

    n_coord = len(xcoords)
    for i in range(n_coord):
        x = xcoords[i]
        y = ycoords[i]
        if (not 0 <= x <= 2048) or (not 0 <= y <= 2048):
            out_vals[i] = 0
        else:
            out_vals[i] = image[x, y]

    return out_vals
