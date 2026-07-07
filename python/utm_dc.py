"""UTM conversion utilities ported from the MATLAB repo.

geog2utm_nodisp : port of geog2utm_nodisp.m (G. Voulgaris, USC 2000)
utm2ll          : port of utm2ll.m (F. Beauducel, IGN algorithms)

Both keep the exact constants and behavior of the MATLAB originals so that
results match the MATLAB eddy_identification_winding repo.
"""

import numpy as np

_UTM_RANGE = np.arange(-180, 181, 6)
_UTM_MERIDIAN = np.arange(-177, 179, 6)


def auto_zone(lon):
    """The UTM zone geog2utm_nodisp picks automatically for these longitudes."""
    lon = np.asarray(lon, dtype=float)
    LL = np.min(lon) + (np.max(lon) - np.min(lon)) / 2
    zone = 0
    for i in range(60):
        if _UTM_RANGE[i] < LL <= _UTM_RANGE[i + 1]:
            zone = i + 1
    return zone


def geog2utm_nodisp(lon, lat, LON=None, LAT=None):
    """Convert lon/lat (deg) to UTM eastings/northings in km.

    Modes (matching the MATLAB nargin behavior):
      geog2utm_nodisp(lon, lat)            -> auto zone from mean lon, 500 km false easting
      geog2utm_nodisp(lon, lat, zone)      -> given integer UTM zone, 500 km false easting
      geog2utm_nodisp(lon, lat, LON, LAT)  -> local origin at (LON, LAT), no false easting

    Returns (X, Y, UTMzone) with X, Y in km as 1-D arrays (or scalars in, arrays out).
    """
    lon = np.atleast_1d(np.asarray(lon, dtype=float)).ravel()
    lat = np.atleast_1d(np.asarray(lat, dtype=float)).ravel()
    UTMzone = 0

    if LAT is None:
        if LON is None:  # auto zone from center longitude
            LL = np.min(lon) + (np.max(lon) - np.min(lon)) / 2
            for i in range(60):
                if _UTM_RANGE[i] < LL <= _UTM_RANGE[i + 1]:
                    UTMzone = i + 1  # 1-based zone number as in MATLAB
        else:
            UTMzone = int(LON)
        XOFF = 500.0
        LAT0 = 0.0
        LON0 = float(_UTM_MERIDIAN[UTMzone - 1])
    else:
        XOFF = 0.0
        LON0 = float(LON)
        LAT0 = float(LAT)

    # constants (Clarke 1866-like, as in the MATLAB original)
    alpha = 6378206.4
    e2 = 0.00676866
    DEG2RAD = 2 * np.pi / 360
    ko = 0.9996

    latr = lat * DEG2RAD
    lonr = lon * DEG2RAD
    LATr = LAT0 * DEG2RAD
    LONr = LON0 * DEG2RAD

    ee2 = e2 / (1 - e2)
    N = alpha / np.sqrt(1 - e2 * np.sin(latr) ** 2)
    T = np.tan(latr) ** 2
    C = ee2 * np.cos(latr) ** 2
    A = (lonr - LONr) * np.cos(latr)
    M = 111132.0894 * latr / DEG2RAD - 16216.94 * np.sin(2 * latr) \
        + 17.21 * np.sin(4 * latr) - 0.02 * np.sin(6 * latr)
    Mo = 111132.0894 * LATr / DEG2RAD - 16216.94 * np.sin(2 * LATr) \
        + 17.21 * np.sin(4 * LATr) - 0.02 * np.sin(6 * LATr)

    X = ko * N * (A + (1 - T + C) * (A ** 3 / 6)
                  + (5 - 18 * T + T ** 2 + 72 * C - 58 * ee2) * (A ** 5 / 120))
    Y = ko * (M - Mo + N * np.tan(latr) * ((A ** 2 / 2)
              + (5 - T + 9 * C + 4 * C ** 2) * (A ** 4 / 24)
              + (61 - 58 * T + T ** 2 + 600 * C - 330 * ee2) * (A ** 6 / 720)))

    X = X / 1000 + XOFF  # km
    Y = Y / 1000         # km
    return X, Y, UTMzone


def _coef(e, m):
    """Projection coefficient polynomials from utm2ll.m."""
    if m == 0:
        c0 = [[-175 / 16384, 0, -5 / 256, 0, -3 / 64, 0, -1 / 4, 0, 1],
              [-105 / 4096, 0, -45 / 1024, 0, -3 / 32, 0, -3 / 8, 0, 0],
              [525 / 16384, 0, 45 / 1024, 0, 15 / 256, 0, 0, 0, 0],
              [-175 / 12288, 0, -35 / 3072, 0, 0, 0, 0, 0, 0],
              [315 / 131072, 0, 0, 0, 0, 0, 0, 0, 0]]
    elif m == 1:
        c0 = [[-175 / 16384, 0, -5 / 256, 0, -3 / 64, 0, -1 / 4, 0, 1],
              [1 / 61440, 0, 7 / 2048, 0, 1 / 48, 0, 1 / 8, 0, 0],
              [559 / 368640, 0, 3 / 1280, 0, 1 / 768, 0, 0, 0, 0],
              [283 / 430080, 0, 17 / 30720, 0, 0, 0, 0, 0, 0],
              [4397 / 41287680, 0, 0, 0, 0, 0, 0, 0, 0]]
    else:
        raise ValueError('m must be 0 or 1')
    return np.array([np.polyval(row, e) for row in c0])


def utm2ll(x, y, zone, datum='wgs84'):
    """Convert UTM x, y in METERS (zone integer, negative for S hemisphere)
    to (lat, lon) in degrees. Port of utm2ll.m (WGS84 default datum)."""
    datums = {'wgs84': (6378137.0, 298.257223563),
              'nad83': (6378137.0, 298.257222101),
              'grs80': (6378137.0, 298.257222101),
              'nad27': (6378206.4, 294.978698214),
              'int24': (6378388.0, 297.0),
              'clk66': (6378206.4, 294.978698214)}
    A1, F1 = datums[datum]

    x = np.atleast_1d(np.asarray(x, dtype=float)).ravel()
    y = np.atleast_1d(np.asarray(y, dtype=float)).ravel()
    f = int(zone)

    D0 = 180 / np.pi
    maxiter = 100
    eps = 1e-11

    K0 = 0.9996
    X0 = 500000.0
    Y0 = 1e7 if f < 0 else 0.0
    P0 = 0.0
    L0 = (6 * abs(f) - 183) / D0
    E1 = np.sqrt((A1 ** 2 - (A1 * (1 - 1 / F1)) ** 2) / A1 ** 2)
    N = K0 * A1

    C = _coef(E1, 0)
    YS = Y0 - N * (C[0] * P0 + C[1] * np.sin(2 * P0) + C[2] * np.sin(4 * P0)
                   + C[3] * np.sin(6 * P0) + C[4] * np.sin(8 * P0))

    C = _coef(E1, 1)
    zt = (y - YS) / N / C[0] + 1j * (x - X0) / N / C[0]
    z = zt - C[1] * np.sin(2 * zt) - C[2] * np.sin(4 * zt) \
        - C[3] * np.sin(6 * zt) - C[4] * np.sin(8 * zt)
    L = np.real(z)
    LS = np.imag(z)

    l = L0 + np.arctan(np.sinh(LS) / np.cos(L))
    p = np.arcsin(np.sin(L) / np.cosh(LS))

    L = np.log(np.tan(np.pi / 4 + p / 2))

    p = 2 * np.arctan(np.exp(L)) - np.pi / 2
    p0 = np.full_like(p, np.nan)
    n = 0
    while np.any(np.isnan(p0) | (np.abs(p - p0) > eps)) and n < maxiter:
        p0 = p
        es = E1 * np.sin(p0)
        p = 2 * np.arctan(((1 + es) / (1 - es)) ** (E1 / 2) * np.exp(L)) - np.pi / 2
        n += 1

    return p * D0, l * D0
