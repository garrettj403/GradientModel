"""Calculate properties of a rough metal surface using the Gradient Model.

Uses the closed-form expression from:

   D. N. Grujic, “Closed-Form Solution of Rough Conductor Surface Impedance,”
   IEEE Trans. Microw. Theory Tech., vol. 66, no. 11, pp. 4677–4683, 2018,
   doi: 10.1109/TMTT.2018.2864586.

"""

import numpy as np
from mpmath import hyp2f1, hyp3f2
from numpy import ndarray, pi, sqrt, exp
from scipy.constants import mu_0


def mag_field(x, f, rq, sigma0=5.8e7):
    """Calculate magnetic field tangential to metal surface.

    Args:
        x: position
        f: frequency
        rq: rms surface roughness
        sigma0: dc conductivity

    Returns:
        magnetic field

    """

    x_is_array = isinstance(x, ndarray) and len(x) > 1
    f_is_array = isinstance(f, ndarray) and len(f) > 1
    if x_is_array and f_is_array:
        print("x and f can't both be arrays")
        raise ValueError

    # Constants
    xi = 0.5
    chi = sqrt(2) * rq

    # Angular frequency
    w = 2 * pi * f
    
    # Eqns 15 and 21
    alpha = (1 + 1j) / 2 * rq * sqrt(mu_0 * w * sigma0)
    beta = 0.5 * (sqrt(1 + 4 * alpha ** 2) - 1)

    # Eqn 32
    zeta = 1 / (1 + exp(2 * (x / chi + xi)))

    # Coefficients
    a1 = alpha + beta
    a2 = alpha - beta - 1
    b1 = 1 + 2 * alpha

    if isinstance(x, ndarray) and len(x) > 1:
        mag = np.empty_like(x, dtype=complex)
        for i, _z in np.ndenumerate(zeta):
            mag[i] = _z ** alpha 
            mag[i] *= hyp2f1(a1, a2, b1, _z)
    elif isinstance(f, ndarray) and len(f) > 1:
        mag = np.empty_like(f, dtype=complex)
        for i in range(len(f)):
            mag[i] = zeta ** alpha[i] * hyp2f1(a1[i], a2[i], b1[i], zeta)
    else:
        mag = zeta ** alpha * hyp2f1(a1, a2, b1, zeta)
    
    return mag


def surface_impedance(f, rq, x0=None, sigma0=5.8e7):
    """Calculate surface impedance of rough metal.

    Args:
        f: frequency
        rq: rms surface roughness
        x0: starting point for integral
        sigma0: dc conductivity

    Returns:
        surface impedance

    """

    f_is_array = isinstance(f, ndarray) and len(f) > 1

    # Constants
    xi = 0.5
    chi = sqrt(2) * rq
    if x0 is None:
        x0 = -5 * rq
        
    # Angular frequency
    w = 2 * pi * f
    
    # Eqns 15 and 21
    # frequency-dependent
    alpha = (1 + 1j) / 2 * rq * sqrt(mu_0 * w * sigma0)
    beta = 0.5 * (sqrt(1 + 4 * alpha ** 2) - 1)

    # Eqn 32
    # position-dependent
    zeta = 1 / (1 + exp(2 * (x0 / chi + xi)))

    # Magnetic field, mag
    mag = mag_field(x0, f, rq, sigma0=sigma0)
    
    # Anti-derivative, bb
    a1 = 1 + alpha - beta
    a2 = 2 + alpha + beta
    a3 = alpha
    b1 = 1 + 2 * alpha
    b2 = 1 + alpha
    if f_is_array:
        f0 = np.empty_like(f, dtype=complex)
        f1 = np.empty_like(f, dtype=complex)
        for i in range(len(f)):
            f1[i] = hyp3f2(a1[i], a2[i], a3[i] + 1, b1[i], b2[i] + 1, zeta)
            f0[i] = hyp3f2(a1[i], a2[i], a3[i], b1[i], b2[i], zeta)
        bb = chi / 2 * (zeta ** alpha) * (zeta / (1 + alpha) * f1 - f0 / alpha)
    else:
        f1 = hyp3f2(a1, a2, a3 + 1, b1, b2 + 1, zeta)
        f0 = hyp3f2(a1, a2, a3, b1, b2, zeta)
        bb = chi / 2 * (zeta ** alpha) * (zeta / (1 + alpha) * f1 - f0 / alpha)
    
    return -1j * mu_0 * w * bb / mag


def rough_properties(f, rq, x0=None, sigma0=5.8e7):
    """Calculate surface properties of rough metal.

    Args:
        f: frequency
        rq: rms surface roughness
        x0: starting point for integral
        sigma0: dc conductivity

    Returns:
        surface impedance, effective conductivity, effective permeability

    """

    # Angular frequency
    w = 2 * pi * f

    # Surface roughness
    zs_rough = surface_impedance(f, rq, x0=x0, sigma0=sigma0)

    # Effective conductivity
    cond_eff = mu_0 * w / (2 * zs_rough.real ** 2)

    # Effective permeability
    mu_eff = 2 * sigma0 * zs_rough.imag ** 2 / w

    return zs_rough, cond_eff, mu_eff / mu_0
