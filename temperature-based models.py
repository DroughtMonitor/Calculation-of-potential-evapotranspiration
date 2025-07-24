__author__ = "Weiqi Liu"
__copyright__ = "Copyright (C) 2025 Weiqi Liu"
__license__ = "NIEER"
__version__ = "2025.07"
__Reference paper__ = "Impact of the potential evapotranspiration models on drought monitoring"
__Description__ =  The temperature-based models.

import numpy as np
import calendar

# Default number of days in each month
_MONTHDAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
_LEAP_MONTHDAYS = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def thornthwaite(monthly_t, monthly_mean_dlh, year=None):
    """
    Estimate monthly PET using Thornthwaite (1948) method.

    Parameters
    ----------
    monthly_t : list or np.ndarray
        Mean daily temperature for each month [°C]
    monthly_mean_dlh : list or np.ndarray
        Mean daily daylight hours for each month [hours]
    year : int, optional
        Year (used for leap year check), default assumes non-leap year

    Returns
    -------
    list
        Monthly potential evapotranspiration [mm/month]
    """
    if len(monthly_t) != 12:
        raise ValueError(f'monthly_t should have 12 values, got {len(monthly_t)}')
    if len(monthly_mean_dlh) != 12:
        raise ValueError(f'monthly_mean_dlh should have 12 values, got {len(monthly_mean_dlh)}')

    month_days = _LEAP_MONTHDAYS if (year and calendar.isleap(year)) else _MONTHDAYS

    adj_t = [max(0, t) for t in monthly_t]  # Negative temps set to 0

    I = sum((t / 5.0) ** 1.514 for t in adj_t if t > 0)
    a = (6.75e-7 * I**3) - (7.71e-5 * I**2) + (1.792e-2 * I) + 0.49239

    pet = []
    for Ta, L, N in zip(adj_t, monthly_mean_dlh, month_days):
        pet_month = 1.6 * (L / 12.0) * (N / 30.0) * ((10.0 * Ta / I) ** a) * 10.0
        pet.append(pet_month)
    return pet


def blaney_criddle(tmean, py, a=-1.55, b=0.96):
    """
    Estimate PET using Blaney-Criddle method.

    Parameters
    ----------
    tmean : float or np.ndarray
        Mean daily air temperature [°C]
    py : float or np.ndarray
        Relative daylength ratio (fraction of max annual)
    a : float
        Calibration coefficient
    b : float
        Calibration coefficient

    Returns
    -------
    float or np.ndarray
        PET [mm/day]
    """
    pet = a + b * (py * (0.457 * tmean + 8.128))
    return clip_zeros(pet)


def kharrufa(tmean, py):
    """
    Estimate PET using Kharrufa method.

    Parameters
    ----------
    tmean : float or np.ndarray
        Mean daily temperature [°C]
    py : float or np.ndarray
        Relative daylength ratio

    Returns
    -------
    float or np.ndarray
        PET [mm/day]
    """
    pet = 0.34 * py * np.power(np.maximum(tmean, 0), 1.3)
    return clip_zeros(pet)


def hamon(tmean, dl):
    """
    Estimate PET using Hamon method.

    Parameters
    ----------
    tmean : float or np.ndarray
        Mean daily temperature [°C]
    dl : float or np.ndarray
        Daylight hours [h]

    Returns
    -------
    float or np.ndarray
        PET [mm/day]
    """
    pet = np.power(dl / 12, 2) * np.exp(tmean / 16)
    return clip_zeros(pet)


def linacre(tmean, elevation, lat, tdew):
    """
    Estimate PET using Linacre method.

    Parameters
    ----------
    tmean : float or np.ndarray
        Mean temperature [°C]
    elevation : float
        Site elevation [m]
    lat : float
        Latitude [°]
    tdew : float or np.ndarray
        Dew point temperature [°C]

    Returns
    -------
    float or np.ndarray
        PET [mm/day]
    """
    th = tmean + 0.006 * elevation
    pet = (500 * th / (100 - lat) + 15 * (tmean - tdew)) / (80 - tmean)
    return clip_zeros(pet)


def romanenko(tmean, rh):
    """
    Estimate PET using Romanenko method.

    Parameters
    ----------
    tmean : float or np.ndarray
        Mean daily temperature [°C]
    rh : float or np.ndarray
        Relative humidity [%]

    Returns
    -------
    float or np.ndarray
        PET [mm/day]
    """
    pet = 4.5 * ((1 + tmean / 25) ** 2) * (1 - rh / 100)
    return clip_zeros(pet)


def schendel(tmean, rh):
    """
    Estimate PET using Schendel method.

    Parameters
    ----------
    tmean : float or np.ndarray
        Mean daily temperature [°C]
    rh : float or np.ndarray
        Relative humidity [%]

    Returns
    -------
    float or np.ndarray
        PET [mm/day]
    """
    pet = 16 * tmean / rh
    return clip_zeros(pet)


def clip_zeros(x):
    """
    Clip negative values to zero.

    Parameters
    ----------
    x : float, list, or np.ndarray

    Returns
    -------
    Same type as input with negatives set to zero.
    """
    return np.maximum(x, 0)
