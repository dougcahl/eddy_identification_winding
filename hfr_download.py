#!/usr/bin/env python3
"""Download HF radar surface currents (HFRNet RTV) from NOAA CoastWatch
ERDDAP and convert them to this toolbox's standard input NetCDF schema
(the same layout as data/data2.nc).

Source: https://coastwatch.pfeg.noaa.gov/erddap/griddap/ (griddap; a plain
HTTP GET returns a NetCDF subset). The near-real-time RTV datasets keep a
ROLLING window of roughly the last 3 months — older requests return an
empty-subset error. The original HFRNet THREDDS server
(hfrnet-tds.ucsd.edu) used by v1.0 is no longer online; an alternative
OPeNDAP mirror exists at
https://dods.ndbc.noaa.gov/thredds/catalog/hfradar/catalog.html.

Usage (standalone):
  python hfr_download.py ucsdHfrW6 2026-07-25 2026-07-27 -125 -120.5 33.5 38.5 \
         -o data/hfr_test.nc [--max-hdop 1.25]
"""

import os
import sys
import time as _time
import tempfile
import urllib.request
import urllib.parse
import urllib.error
from datetime import datetime, timezone

import numpy as np
from netCDF4 import Dataset, num2date

ERDDAP_BASE = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap'

DATASETS = {  # dataset id -> nominal resolution (km)
    'ucsdHfrE1': 1, 'ucsdHfrE2': 2, 'ucsdHfrE6': 6,
    'ucsdHfrW1': 1, 'ucsdHfrW2': 2, 'ucsdHfrW6': 6,
}
REGION_PRESETS = {'USEGC': 'E', 'USWC': 'W'}

DATENUM_1970 = 719529.0  # MATLAB datenum of 1970-01-01
ERDDAP_MAX_BYTES = 1.5e9  # stay under ERDDAP's ~2 GB response cap


def resolve_dataset(dataset, res_km=None):
    """Accept a dataset id ('ucsdHfrW6') or a region preset ('USWC') plus
    res_km (1|2|6); return (dataset_id, res_km)."""
    if dataset in DATASETS:
        return dataset, DATASETS[dataset]
    if dataset in REGION_PRESETS:
        if res_km not in (1, 2, 6):
            raise ValueError('preset %r needs res_km in (1, 2, 6)' % dataset)
        did = 'ucsdHfr%s%d' % (REGION_PRESETS[dataset], res_km)
        return did, res_km
    raise ValueError('unknown dataset %r; valid ids: %s; presets: %s'
                     % (dataset, sorted(DATASETS), sorted(REGION_PRESETS)))


def _iso(t):
    if isinstance(t, datetime):
        return t.strftime('%Y-%m-%dT%H:%M:%SZ')
    t = str(t)
    for fmt in ('%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M', '%Y-%m-%d'):
        try:
            return datetime.strptime(t.rstrip('Z'), fmt).strftime('%Y-%m-%dT%H:%M:%SZ')
        except ValueError:
            continue
    raise ValueError('cannot parse time %r (use YYYY-MM-DD[THH:MM[:SS]])' % t)


def _build_query(dataset_id, t0_iso, t1_iso, bbox, want_hdop):
    lon_min, lon_max, lat_min, lat_max = bbox
    sel = '[(%s):(%s)][(%.6g):(%.6g)][(%.6g):(%.6g)]' % (
        t0_iso, t1_iso, lat_min, lat_max, lon_min, lon_max)
    variables = ['water_u', 'water_v'] + (['hdop'] if want_hdop else [])
    return ','.join(v + sel for v in variables)


def _fetch(url, tmp_path, timeout, retries):
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp, \
                    open(tmp_path, 'wb') as f:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return
        except urllib.error.HTTPError as e:
            body = e.read().decode(errors='replace')
            if 'nRows = 0' in body or 'no matching results' in body.lower():
                raise RuntimeError(
                    'ERDDAP returned an empty subset. The RTV near-real-time '
                    'datasets keep a rolling ~3-month window; check the '
                    'current coverage at %s.html' % url.split('.nc?')[0]) from None
            if e.code == 404 and 'Resource not found' in body:
                raise ValueError('unknown ERDDAP dataset; valid ids: %s'
                                 % sorted(DATASETS)) from None
            raise RuntimeError('ERDDAP error HTTP %d: %s'
                               % (e.code, body[:500])) from None
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
            if attempt < retries - 1:
                _time.sleep(2 ** attempt * 5)
    raise RuntimeError('download failed after %d attempts: %s' % (retries, last))


def _axis(nc, *names):
    for n in names:
        if n in nc.variables:
            return n
    raise KeyError('none of %s found in downloaded file (has: %s)'
                   % (names, list(nc.variables)))


def _convert(tmp_path, out_nc, dataset_id, url, max_hdop):
    with Dataset(tmp_path) as src:
        latn = _axis(src, 'latitude', 'lat')
        lonn = _axis(src, 'longitude', 'lon')
        un = _axis(src, 'water_u', 'u')
        vn = _axis(src, 'water_v', 'v')
        lat = np.array(src.variables[latn][:], dtype=float)
        lon = np.array(src.variables[lonn][:], dtype=float)
        tvar = src.variables['time']
        tvals = np.array(tvar[:], dtype=float)
        # robust epoch conversion via the units attribute
        dts = num2date(tvals, tvar.units,
                       only_use_cftime_datetimes=False,
                       only_use_python_datetimes=True)
        time_days = np.array(
            [(d.replace(tzinfo=None) - datetime(1970, 1, 1)).total_seconds()
             for d in dts]) / 86400.0

        u = np.ma.filled(src.variables[un][:], np.nan).astype(float)
        v = np.ma.filled(src.variables[vn][:], np.nan).astype(float)
        hdop = None
        if max_hdop is not None:
            hn = _axis(src, 'hdop')
            hdop = np.ma.filled(src.variables[hn][:], np.nan).astype(float)

    # normalize longitudes to -180..180 and force ascending axes
    lon = ((lon + 180) % 360) - 180
    if lon[0] > lon[-1]:
        lon = lon[::-1]
        u = u[:, :, ::-1]
        v = v[:, :, ::-1]
        if hdop is not None:
            hdop = hdop[:, :, ::-1]
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        u = u[:, ::-1, :]
        v = v[:, ::-1, :]
        if hdop is not None:
            hdop = hdop[:, ::-1, :]

    # sort/dedupe time
    order = np.argsort(time_days)
    time_days = time_days[order]
    u = u[order]
    v = v[order]
    if hdop is not None:
        hdop = hdop[order]
    keep = np.concatenate([[True], np.diff(time_days) > 1e-9])
    time_days = time_days[keep]
    u = u[keep]
    v = v[keep]
    if hdop is not None:
        hdop = hdop[keep]

    if max_hdop is not None:
        bad = ~(hdop <= max_hdop)  # NaN hdop also fails QC
        u[bad] = np.nan
        v[bad] = np.nan

    valid_per_step = np.isfinite(u).sum(axis=(1, 2))
    if valid_per_step.max() == 0:
        raise RuntimeError('every timestep is all-NaN: the bounding box has '
                           'no radar coverage in this window')
    n_empty = int((valid_per_step == 0).sum())
    if n_empty:
        print('warning: %d of %d timesteps have no valid vectors'
              % (n_empty, len(time_days)))

    # (time, lat, lon) -> (time, i=lon, j=lat), matching data2.nc
    u = np.swapaxes(u, 1, 2)
    v = np.swapaxes(v, 1, 2)
    lats, lons = np.meshgrid(lat, lon)  # (nlon, nlat)

    with Dataset(out_nc, 'w', format='NETCDF4') as nc:
        nc.title = 'HF radar surface currents (HFRNet RTV via ERDDAP)'
        nc.source = dataset_id
        nc.request_url = url
        nc.history = 'downloaded ' + datetime.now(timezone.utc).strftime(
            '%Y-%m-%dT%H:%M:%SZ')
        if max_hdop is not None:
            nc.max_hdop = float(max_hdop)

        nc.createDimension('i', lons.shape[0])
        nc.createDimension('j', lons.shape[1])
        nc.createDimension('time', len(time_days))

        vlon = nc.createVariable('lon', 'f8', ('i', 'j'), zlib=True)
        vlon.units = 'degrees_east'
        vlat = nc.createVariable('lat', 'f8', ('i', 'j'), zlib=True)
        vlat.units = 'degrees_north'
        vt = nc.createVariable('time', 'f8', ('time',))
        vt.units = 'days since 1970-01-01 00:00:00'
        vt.calendar = 'standard'
        vdn = nc.createVariable('matlab_datenum', 'f8', ('time',))
        vdn.long_name = 'MATLAB datenum (days since year 0)'
        vu = nc.createVariable('u', 'f8', ('time', 'i', 'j'), zlib=True,
                               fill_value=np.nan)
        vu.units = 'm s-1'
        vu.long_name = 'eastward surface current velocity'
        vv = nc.createVariable('v', 'f8', ('time', 'i', 'j'), zlib=True,
                               fill_value=np.nan)
        vv.units = 'm s-1'
        vv.long_name = 'northward surface current velocity'

        vlon[:] = lons
        vlat[:] = lats
        vt[:] = time_days
        vdn[:] = time_days + DATENUM_1970
        vu[:] = u
        vv[:] = v

    print('wrote %s  (%d x %d grid, %d timesteps, %.0f%% mean coverage)'
          % (out_nc, lons.shape[0], lons.shape[1], len(time_days),
             100 * np.isfinite(u).mean()))


def download_hfr(dataset, t_start, t_end, bbox, out_nc, max_hdop=None,
                 timeout=120, retries=3, base_url=ERDDAP_BASE, res_km=None):
    """Download water_u/water_v (and hdop when max_hdop is set) for
    bbox=(lon_min, lon_max, lat_min, lat_max) between t_start and t_end and
    write the toolbox-standard input NetCDF out_nc. Returns out_nc."""
    dataset_id, res_km = resolve_dataset(dataset, res_km)

    lon_min, lon_max, lat_min, lat_max = [float(x) for x in bbox]
    lon_min = ((lon_min + 180) % 360) - 180
    lon_max = ((lon_max + 180) % 360) - 180
    if lon_min >= lon_max:
        raise ValueError('bbox crosses the dateline or lon_min >= lon_max '
                         'after normalizing to -180..180 (unsupported)')
    if lat_min >= lat_max:
        raise ValueError('lat_min must be < lat_max')
    bbox = (lon_min, lon_max, lat_min, lat_max)

    t0, t1 = _iso(t_start), _iso(t_end)
    if t0 >= t1:
        raise ValueError('t_start must be before t_end')

    # rough size estimate against the ERDDAP response cap
    step_deg = res_km / 111.0
    n_hours = (datetime.strptime(t1, '%Y-%m-%dT%H:%M:%SZ')
               - datetime.strptime(t0, '%Y-%m-%dT%H:%M:%SZ')).total_seconds() / 3600
    est = ((lon_max - lon_min) / step_deg) * ((lat_max - lat_min) / step_deg) \
        * max(n_hours, 1) * 8 * (3 if max_hdop is not None else 2)
    if est > ERDDAP_MAX_BYTES:
        raise ValueError('request too large (~%.1f GB estimated); split the '
                         'time range into smaller chunks' % (est / 1e9))

    query = _build_query(dataset_id, t0, t1, bbox, max_hdop is not None)
    url = '%s/%s.nc?%s' % (base_url, dataset_id,
                           urllib.parse.quote(query, safe='[]():,-'))

    os.makedirs(os.path.dirname(os.path.abspath(out_nc)), exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(suffix='.nc',
                                    dir=os.path.dirname(os.path.abspath(out_nc)))
    os.close(fd)
    try:
        print('downloading %s %s to %s ...' % (dataset_id, t0, t1))
        _fetch(url, tmp_path, timeout, retries)
        with open(tmp_path, 'rb') as f:
            magic = f.read(4)
        if not (magic.startswith(b'CDF') or magic.startswith(b'\x89HDF')):
            with open(tmp_path, 'r', errors='replace') as f:
                raise RuntimeError('server returned a non-NetCDF response: '
                                   + f.read(500))
        _convert(tmp_path, out_nc, dataset_id, url, max_hdop)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return out_nc


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('dataset', help='ERDDAP id (e.g. ucsdHfrW6) or preset '
                                   '(USEGC/USWC, with --res)')
    p.add_argument('start')
    p.add_argument('end')
    p.add_argument('bbox', nargs=4, type=float,
                   metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'))
    p.add_argument('-o', '--out', default=None)
    p.add_argument('--res', type=int, default=None, choices=(1, 2, 6))
    p.add_argument('--max-hdop', type=float, default=None)
    a = p.parse_args()
    did, _ = resolve_dataset(a.dataset, a.res)
    out = a.out or 'data/%s_%s.nc' % (did, a.start.replace('-', ''))
    download_hfr(a.dataset, a.start, a.end, a.bbox, out,
                 max_hdop=a.max_hdop, res_km=a.res)
