"""Triaxial (ellipsoidal) projected halo profiles and associated functions.
"""

import cluster_toolkit
from cluster_toolkit import _dcast
import numpy as np

def mapped_radii(z, R, i, phi, q, s):
    r"""Radial vector norm in the frame of the triaxial
    halo, mapped to the isodensity contour of a spherical halo
    at r = \sqrt(r^2 + z^2). Units are Mpc/h comoving for distances
    and radians for angles.

    Args:
        z (float or array like): vertical distance above the x-y plane
        R (float or array like): projected distance in the x-y plane
        i (float or array like): orientation angle in radians
        phi (float or array like): azimuthal angle
        q (float): minor-to-major axis ratio
        s (float): intermediate-to-major axis ratio

    Returns:
       float or array like: mapped radii

    """
    if q < 0 or q > 1:
        raise Exception("q must be > 0 and <= 1.")

    if s < q or s > 1:
        raise Exception("s must be >= q and <= 1.")

    z = np.asarray(z)
    R = np.asarray(R)
    i = np.asarray(i)
    phi = np.asarray(phi)
    squeezed_dims = []
    if z.ndim == 0:
        z = z[None]
        squeezed_dims.append(0)
    if R.ndim == 0:
        R = R[None]
        squeezed_dims.append(1)
    if i.ndim == 0:
        i = i[None]
        squeezed_dims.append(2)
    if phi.ndim == 0:
        phi = phi[None]
        squeezed_dims.append(3)

    if z.ndim >1 or R.ndim > 1 or i.ndim > 1 or phi.ndim > 1:
        raise Exception("Position inputs cannot be more than "+
                        "one dimensional arrays.")

    mapped_r = np.zeros((len(z), len(R), len(i), len(phi)))
    for j in range(len(z)):
        for k in range(len(R)):
            for l in range(len(i)):
                sin_i = np.sin(i[l])
                cos_i = np.cos(i[l])
                for m in range(len(phi)):
                    sin_phi = np.sin(phi[m])
                    cos_phi = np.cos(phi[m])
                    mapped_r[j,k,l,m] =\
                        cluster_toolkit._lib.mapped_radii(z[j], R[k],
                                                          sin_i, cos_i,
                                                          sin_phi, cos_phi,
                                                          q, s)

    #If we have to squeeze any dimensions, do so and return
    if squeezed_dims:
        return np.squeeze(mapped_r, tuple(squeezed_dims))
    return mapped_r

