__author__ = "Weiqi Liu"
__copyright__ = "Copyright (C) 2025 Weiqi Liu"
__license__ = "NIEER"
__version__ = "2025.07"
__Reference paper__ = "Impact of the potential evapotranspiration models on drought monitoring"
__Description__ =  The Combination models.

import numpy as np

def fao56_penman_monteith(net_rad: float, t: float, ws: float,
                           svp: float, avp: float, delta_svp: float,
                           psy: float, shf: float = 0.0) -> float:
    """
    FAO-56 Penman-Monteith equation for reference evapotranspiration (ETo).

    Parameters
    ----------
    net_rad : float
        Net radiation at crop surface [MJ m⁻² day⁻¹]
    t : float
        Air temperature at 2 m height [°C]
    ws : float
        Wind speed at 2 m height [m s⁻¹]
    svp : float
        Saturation vapour pressure [kPa]
    avp : float
        Actual vapour pressure [kPa]
    delta_svp : float
        Slope of saturation vapour pressure curve [kPa °C⁻¹]
    psy : float
        Psychrometric constant [kPa °C⁻¹]
    shf : float, optional
        Soil heat flux [MJ m⁻² day⁻¹], by default 0.0

    Returns
    -------
    float
        Reference evapotranspiration [mm day⁻¹]
    """
    term1 = 0.408 * (net_rad - shf) * delta_svp / (delta_svp + psy * (1 + 0.34 * ws))
    term2 = 900 * ws * (svp - avp) * psy / ((t + 273) * (delta_svp + psy * (1 + 0.34 * ws)))
    pet = term1 + term2
    return clip_zeros(pet)


def penman(lambd: float, net_rad: float, wind: float, delta_svp: float,
           svp: float, avp: float, psy: float,
           aw: float = 1, bw: float = 0.537, ku: float = 6.43,
           shf: float = 0.0) -> float:
    """
    Standard Penman equation for potential evapotranspiration.

    Returns PET [mm day⁻¹]
    """
    fu = ku * (aw + bw * wind)
    denominator = lambd * (delta_svp + psy)
    term1 = delta_svp * (net_rad - shf) / denominator
    term2 = psy * (svp - avp) * fu / denominator
    pet = term1 + term2
    return clip_zeros(pet)


def FAO24_penman(lambd: float, net_rad: float, wind: float, delta_svp: float,
                 svp: float, avp: float, psy: float,
                 aw: float = 1, bw: float = 0.864, ku: float = 2.7,
                 shf: float = 0.0) -> float:
    """
    FAO-24 version of Penman equation.

    Returns PET [mm day⁻¹]
    """
    fu = ku * (aw + bw * wind)
    denominator = lambd * (delta_svp + psy)
    term1 = delta_svp * (net_rad - shf) / denominator
    term2 = psy * (svp - avp) * fu / denominator
    pet = term1 + term2
    return clip_zeros(pet)


def FAO_ppp17_Penman(lambd: float, net_rad: float, wind: float,
                     delta_svp: float, svp: float, avp: float,
                     psy: float, Tmax: float, Tmin: float,
                     aw: float = 1, ku: float = 6.43, shf: float = 0.0) -> float:
    """
    FAO PPP17 version of Penman with temperature-dependent wind coefficient.

    Returns PET [mm day⁻¹]
    """
    delta_T = Tmax - Tmin
    if delta_T < 12:
        bw = 0.54
    else:
        bw = 0.54 + 0.35 * (delta_T - 12) / 4

    fu = ku * (aw + bw * wind)
    denominator = lambd * (delta_svp + psy)
    term1 = delta_svp * (net_rad - shf) / denominator
    term2 = psy * (svp - avp) * fu / denominator
    pet = term1 + term2
    return clip_zeros(pet)


def kimberly_penman(ws: float, net_rad: float, jj: int, avp: float,
                    svp: float, lambd: float, delta_svp: float,
                    psy: float, g: float = 0.0, ku: float = 2.62) -> float:
    """
    Kimberly Penman equation.

    Parameters
    ----------
    jj : int
        Day of year (DOY)
    
    Returns PET [mm day⁻¹]
    """
    # Empirical seasonal wind function
    w = ((0.4 + 1.4 * np.exp(-((jj - 173) / 58) ** 2)) +
         (0.605 + 0.345 * np.exp(-((jj - 243) / 80) ** 2))) * ws

    denominator = lambd * (delta_svp + psy)
    term1 = delta_svp * (net_rad - g) / denominator
    term2 = ku * psy * (svp - avp) * w / denominator
    pet = term1 + term2
    return clip_zeros(pet)
