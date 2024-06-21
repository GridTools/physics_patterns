import numpy as np


def vertical_loop_pattern(
    zsoil,
    nroot,
    sh2o,
    smcwlt,
    smcref,
):
    '''
    This is a pattern in the land surface model where a calculation is extended
    below the surface by a certain number of layers based on something like
    root depth. The number of vertical levels is different for different grid
    points, but is known at compile time. A similar motif exists in the PBL
    scheme where the calculation is to within the PBL height, but that can
    change at runtime
    '''
    rcsoil = np.zeros(nroot.shape)
    for i in range(zsoil.shape[0]):
        for j in range(zsoil.shape[1]):
            for k in range(nroot[i, j]):
                gx = (sh2o[i, j, k] - smcwlt[i, j]) / (
                    smcref[i, j] - smcwlt[i, j]
                )
                gx = max(0.0, min(1.0, gx))
                if k == 0:
                    rcsoil[i, j] += (zsoil[i, j, 0]/zsoil[i, j, nroot]) * gx
                else:
                    rcsoil[i, j] += (
                        (zsoil[k] - zsoil(k-1)) / zsoil(nroot)
                    ) * gx
    return rcsoil
