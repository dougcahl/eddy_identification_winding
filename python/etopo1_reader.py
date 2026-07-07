"""Read a lon/lat subregion from the ETOPO1 bathymetry binary
(etopo1_ice_g_i2.bin, grid-registered, 1 arc-minute, int16 little-endian) —
the same file m_map's m_etopo2 uses in the MATLAB repo.

Grid layout (from m_etopo2.m): 10801 rows x 21601 cols, row 0 = 90N,
col 0 = 180W, 60 points per degree.
"""

import os
import numpy as np

PTSPERDEG = 60
NCOLS = 360 * PTSPERDEG + 1   # 21601 (grid registration)
NROWS = 180 * PTSPERDEG + 1   # 10801

_HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH_PATHS = [
    os.path.join(_HERE, 'etopo1', 'etopo1_ice_g_i2.bin'),
    os.path.join(_HERE, '..', 'm_map', 'etopo1', 'etopo1_ice_g_i2.bin'),  # repo layout
    os.path.join(_HERE, '..', 'etopo1', 'etopo1_ice_g_i2.bin'),
    r'D:\sync\MATLAB\m_map\etopo1\etopo1_ice_g_i2.bin',
]


def read_etopo1(lon_min, lon_max, lat_min, lat_max, path=None):
    """Return (elev, lons, lats) for the requested box.
    elev is a 2-D int16 array (nlat, nlon) with lats ascending."""
    if path is None:
        for p in SEARCH_PATHS:
            if os.path.isfile(p):
                path = p
                break
    if path is None or not os.path.isfile(path):
        raise FileNotFoundError('etopo1_ice_g_i2.bin not found; put it in ./etopo1/ '
                                'or pass path=')

    data = np.memmap(path, dtype='<i2', mode='r', shape=(NROWS, NCOLS))

    j0 = int(np.floor((lon_min + 180) * PTSPERDEG))
    j1 = int(np.ceil((lon_max + 180) * PTSPERDEG))
    i0 = int(np.floor((90 - lat_max) * PTSPERDEG))
    i1 = int(np.ceil((90 - lat_min) * PTSPERDEG))
    j0 = max(j0, 0); j1 = min(j1, NCOLS - 1)
    i0 = max(i0, 0); i1 = min(i1, NROWS - 1)

    elev = np.array(data[i0:i1 + 1, j0:j1 + 1])
    lons = j0 / PTSPERDEG - 180 + np.arange(elev.shape[1]) / PTSPERDEG
    lats = 90 - i0 / PTSPERDEG - np.arange(elev.shape[0]) / PTSPERDEG

    # flip so latitude ascends
    return elev[::-1, :], lons, lats[::-1]
