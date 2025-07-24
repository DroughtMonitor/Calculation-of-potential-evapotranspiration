__author__ = "Weiqi Liu"
__copyright__ = "Copyright (C) 2025 Weiqi Liu"
__license__ = "NIEER"
__version__ = "2025.07"
__Reference paper__ = "Impact of the potential evapotranspiration models on drought monitoring"
__Description__ =  The radiation-based models.

# This is a code refinement pass for the user's PET models. We will clean up:
# - Typo in docstrings (e.g. lambd：-> lambda:, "Returns：" -> "Returns")
# - Use of print statements
# - Consistent parameter naming and formatting
# - Comment structure
# - Ensure code is ready for GitHub publishing

import numpy as np

def clip_zeros(s):
    """Replace negative values with 0 for numpy arrays or scalars."""
    return np.where(s < 0, 0, s)

def priestley_taylor(rn, dlt, gamma, lambd, alpha=1.26, g=0):
    """Priestley-Taylor method for PET estimation.

    Parameters
    ----------
    rn : float or array
        Net radiation [MJ m-2 d-1]
    dlt : float or array
        Slope of saturation vapour pressure curve [kPa °C⁻¹]
    gamma : float or array
        Psychrometric constant [kPa °C⁻¹]
    lambd : float or array
        Latent heat of vaporization [MJ kg⁻¹]
    alpha : float, optional
        Priestley-Taylor coefficient [-]
    g : float or array, optional
        Soil heat flux [MJ m-2 d-1]

    Returns
    -------
    pet : float or array
        Potential evapotranspiration [mm d⁻¹]
    """
    pet = (alpha * dlt * (rn - g)) / (lambd * (dlt + gamma))
    return clip_zeros(pet)

def makkink(rs, dlt, gamma, lambd, alpha=0.61, b=-0.012):
    """Makkink method for PET estimation."""
    pet = alpha * dlt * rs / ((dlt + gamma) * lambd) + b
    return clip_zeros(pet)

def jensen_haise(tmean, rs, lambd, cr=0.025, tx=-3):
    """Jensen-Haise method for PET estimation."""
    pet = rs / lambd * cr * (tmean - tx)
    return clip_zeros(pet)

def turc(tmean, rs, rh, k=0.013):
    """Turc method for PET estimation (converted radiation units)."""
    rs = rs * 23.885  # MJ/m2/day -> cal/cm2/day
    c = np.where(rh >= 50, 1, 1 + (50 - rh) / 70)
    pet = k * c * tmean * (rs + 50) / (tmean + 15)
    pet = np.where(tmean <= 0, 0, pet)
    return clip_zeros(pet)

def mcguinness_bordne(tmean, ra, lambd, k=0.0147):
    """McGuinness-Bordne method for PET estimation."""
    pet = k * ra * (tmean + 5) / lambd
    return clip_zeros(pet)

def abtew(tmax, rs, lambd, k=0.01786):
    """Abtew method for PET estimation."""
    pet = k * rs * tmax / lambd
    return clip_zeros(pet)

def doorenbos_pruitt(rs, delta, rh, wind, psy, lambd):
    """Doorenbos-Pruitt method for PET estimation."""
    b = -0.3
    a = (1.066 - 0.0013 * rh - 0.0002 * rh * wind +
         0.045 * wind - 0.0000315 * rh ** 2 - 0.0011 * wind ** 2)
    pet = a * (delta / (delta + psy)) * rs / lambd + b
    return clip_zeros(pet)

def irmak(rs, tmean):
    """Irmak method for PET estimation."""
    pet = -0.611 + 0.149 * rs + 0.079 * tmean
    return clip_zeros(pet)

def oudin(tmean, ra, lambd, k1=100, k2=5):
    """Oudin method for PET estimation."""
    pet = ra * (tmean + k2) / lambd / k1
    pet = np.where((tmean + k2) > 0, pet, 0)
    return clip_zeros(pet)

def hargreaves(tmin, tmax, tmean, et_rad):
    """Hargreaves method for PET estimation."""
    pet = 0.0023 * (tmean + 17.8) * np.sqrt(tmax - tmin) * 0.408 * et_rad
    return clip_zeros(pet)

def droogers_allen(tmin, tmax, tmean, et_rad):
    """Droogers-Allen variation of Hargreaves method."""
    pet = 0.0025 * (tmean + 16.8) * np.sqrt(tmax - tmin) * 0.408 * et_rad
    return clip_zeros(pet)

def allen(tmin, tmax, tmean, et_rad):
    """Allen variation of Hargreaves method."""
    pet = 0.003 * (tmean + 20) * ((tmax - tmin) ** 0.4) * 0.408 * et_rad
    return clip_zeros(pet)

def dorji(tmin, tmax, tmean, et_rad):
    """Dorji method for PET estimation."""
    pet = 0.002 * (tmean + 33.9) * ((tmax - tmin) ** 0.296) * 0.408 * et_rad
    return clip_zeros(pet)

