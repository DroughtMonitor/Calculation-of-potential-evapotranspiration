__author__ = "Weiqi Liu"
__copyright__ = "Copyright (C) 2025 Weiqi Liu"
__license__ = "NIEER"
__version__ = "2025.07"
__Reference paper__ = "Impact of the potential evapotranspiration models on drought monitoring"
__Description__ =  The mass-transfer-based models.

def Dalton(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Dalton equation.

    Parameters
    ----------
    wind : float or array-like
        Wind speed at 2 m height [m s-1].
    es : float or array-like
        Saturation vapour pressure [kPa].
    ea : float or array-like
        Actual vapour pressure [kPa].

    Returns
    -------
    pet : float or array-like
        Reference evapotranspiration [mm day-1].
    """
    pet = (0.3648 + 0.07223 * wind) * (es - ea)
    return clip_zeros(pet)


def Trabert(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Trabert equation.
    """
    pet = 0.3075 * (wind ** 0.5) * (es - ea)
    return clip_zeros(pet)


def Meyer(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Meyer equation.
    """
    pet = (0.375 + 0.05026 * wind) * (es - ea)
    return clip_zeros(pet)


def Rohwer(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Rohwer equation.
    """
    pet = 0.44 * (1 + 0.27 * wind) * (es - ea)
    return clip_zeros(pet)


def Penman_mass(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Penman mass transfer method.
    """
    pet = 0.35 * (1 + 0.98 / (100 * wind)) * (es - ea)
    return clip_zeros(pet)


def Albrecht(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Albrecht equation.
    """
    pet = (0.1005 + 0.297 * wind) * (es - ea)
    return clip_zeros(pet)


def Brockamp_Wenner(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Brockamp-Wenner equation.
    """
    pet = 0.543 * (wind ** 0.456) * (es - ea)
    return clip_zeros(pet)


def WMO(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the WMO equation.
    """
    pet = (0.1298 + 0.0934 * wind) * (es - ea)
    return clip_zeros(pet)


def Mahringer(wind, es, ea):
    """
    Estimate reference evapotranspiration (ETo) using the Mahringer equation.
    """
    pet = 0.15072 * ((3.6 * wind) ** 0.5) * (es - ea)
    return clip_zeros(pet)
