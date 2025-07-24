__author__ = "Weiqi Liu"
__copyright__ = "Copyright (C) 2025 Weiqi Liu"
__license__ = "NIEER"
__version__ = "2025.07"
__Reference paper__ = "Impact of the potential evapotranspiration models on drought monitoring","Allen et al., FAO Irrigation and Drainage Paper No. 56 (FAO56)"
__Description__ =  A collection of functions for estimating actual and saturation vapour pressure,
             psychrometric constant, atmospheric pressure, and related parameters as per 
             FAO-56 guidelines.

import math
import numpy as np


def avp_from_tdew(tdew: float) -> float:
    """
    Estimate actual vapour pressure (ea) from dew point temperature.

    Parameters
    ----------
    tdew : float
        Dew point temperature [°C]

    Returns
    -------
    float
        Actual vapour pressure [kPa]
    """
    return 0.6108 * math.exp((17.27 * tdew) / (tdew + 237.3))


def svp_from_t(t: float) -> float:
    """
    Estimate saturation vapour pressure (es) from air temperature.

    Parameters
    ----------
    t : float
        Air temperature [°C]

    Returns
    -------
    float
        Saturation vapour pressure [kPa]
    """
    return 0.6108 * math.exp((17.27 * t) / (t + 237.3))


def delta_svp(t: float) -> float:
    """
    Estimate the slope of the saturation vapour pressure curve at a given temperature.

    Parameters
    ----------
    t : float
        Air temperature [°C]

    Returns
    -------
    float
        Slope of saturation vapour pressure curve [kPa °C⁻¹]
    """
    es = svp_from_t(t)
    return (4098 * es) / ((t + 237.3) ** 2)


def psy_const(atmos_pres: float) -> float:
    """
    Calculate the psychrometric constant.

    Parameters
    ----------
    atmos_pres : float
        Atmospheric pressure [kPa]

    Returns
    -------
    float
        Psychrometric constant [kPa °C⁻¹]
    """
    return 0.000665 * atmos_pres


def atm_pressure(altitude: float) -> float:
    """
    Estimate atmospheric pressure from altitude.

    Parameters
    ----------
    altitude : float
        Elevation above sea level [m]

    Returns
    -------
    float
        Atmospheric pressure [kPa]
    """
    tmp = (293.0 - (0.0065 * altitude)) / 293.0
    return (tmp ** 5.26) * 101.3


def calc_lambda(tmean: float) -> float:
    """
    Calculate latent heat of vaporization (lambda).

    Parameters
    ----------
    tmean : float
        Mean daily air temperature [°C]

    Returns
    -------
    float
        Latent heat of vaporization [MJ kg⁻¹]
    """
    return 2.501 - 0.002361 * tmean


def wind_speed_2m(ws: float, z: float) -> float:
    """
    Convert wind speed measured at height z to wind speed at 2 meters above ground level.

    Parameters
    ----------
    ws : float
        Wind speed measured at height z [m s⁻¹]
    z : float
        Height of wind measurement above the ground [m]

    Returns
    -------
    float
        Wind speed adjusted to 2 meters height [m s⁻¹]
    """
    return ws * (4.87 / math.log((67.8 * z) - 5.42))


def clip_zeros(arr: np.ndarray) -> np.ndarray:
    """
    Replace negative values with 0.

    Parameters
    ----------
    arr : np.ndarray
        Input array (e.g., from Pandas Series or xarray DataArray)

    Returns
    -------
    np.ndarray
        Array with negative values clipped to 0
    """
    return np.where(arr < 0, 0, arr)
