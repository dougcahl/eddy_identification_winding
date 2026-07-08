"""Native NetCDF I/O for eddy identification results and tracks.

The pipeline reads and writes only .nc files. Schema matches convert_data.py:
scalar-per-eddy variables plus ragged streamline/track series stored as
concatenated arrays with *_start/*_len index vectors. Time is CF-style
days since 1970-01-01, with the original MATLAB datenum in `matlab_datenum`.
"""

import os

import numpy as np
from netCDF4 import Dataset

DATENUM_1970 = 719529.0  # MATLAB datenum of 1970-01-01

IDENT_KEYS = ['eddy_center_lat', 'eddy_center_lon', 'eddy_dir',
              'eddy_angular_vel', 'eddy_streamlines', 'eddy_length_x',
              'eddy_length_y', 'eddy_ellipse_theta']

_IDENT_META = {
    'eddy_center_lat': ('degrees_north', ''),
    'eddy_center_lon': ('degrees_east', ''),
    'eddy_dir': ('1', '+1 CCW/cyclonic, -1 CW/anticyclonic'),
    'eddy_angular_vel': ('degree s-1', ''),
    'eddy_streamlines': ('1', 'number of streamlines in the eddy cluster'),
    'eddy_length_x': ('km', 'ellipse semi-axis 1'),
    'eddy_length_y': ('km', 'ellipse semi-axis 2'),
    'eddy_ellipse_theta': ('degrees', 'ellipse rotation, math notation'),
}


def save_ident_nc(path, time_datenum, UTMzone, scalars, stream_lat, stream_lon,
                  stream_eddy):
    """Write one timestep's identification results.

    scalars     : dict of the IDENT_KEYS -> 1-D arrays (len = n eddies)
    stream_lat/lon : list of 1-D arrays, one per streamline
    stream_eddy : 0-based eddy index for each streamline
    """
    if not path.endswith('.nc'):
        path = path + '.nc'
    n = len(np.ravel(scalars['eddy_streamlines']))
    with Dataset(path, 'w', format='NETCDF4') as nc:
        nc.title = 'winding-angle eddy identification results'
        nc.UTMzone = int(UTMzone)
        nc.createDimension('eddy', n)
        vt = nc.createVariable('time', 'f8', ())
        vt.units = 'days since 1970-01-01 00:00:00'
        vt[...] = float(time_datenum) - DATENUM_1970
        vdn = nc.createVariable('matlab_datenum', 'f8', ())
        vdn[...] = float(time_datenum)

        for k in IDENT_KEYS:
            units, ln = _IDENT_META[k]
            v = nc.createVariable(k, 'f8', ('eddy',))
            v.units = units
            if ln:
                v.long_name = ln
            v[:] = np.ravel(scalars[k]).astype(float) if n else []

        nstreams = len(stream_lat)
        npts = int(sum(len(s) for s in stream_lat))
        nc.createDimension('stream', nstreams)
        nc.createDimension('stream_pt', npts)
        v1 = nc.createVariable('stream_lat', 'f8', ('stream_pt',), zlib=True)
        v2 = nc.createVariable('stream_lon', 'f8', ('stream_pt',), zlib=True)
        v3 = nc.createVariable('stream_eddy', 'i4', ('stream',))
        v3.long_name = '0-based eddy index each streamline belongs to'
        v4 = nc.createVariable('stream_start', 'i4', ('stream',))
        v4.long_name = '0-based start offset of each streamline in stream_lat/lon'
        v5 = nc.createVariable('stream_len', 'i4', ('stream',))
        if nstreams:
            v1[:] = np.concatenate([np.ravel(s) for s in stream_lat])
            v2[:] = np.concatenate([np.ravel(s) for s in stream_lon])
            v3[:] = np.asarray(stream_eddy, dtype=np.int32)
            lens = np.array([len(s) for s in stream_lat], dtype=np.int32)
            v4[:] = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int32)
            v5[:] = lens


def load_ident_nc(path):
    """Read one timestep's results -> dict with IDENT_KEYS arrays, 'time'
    (MATLAB datenum), 'UTMzone', and streamlines as lists of arrays."""
    with Dataset(path) as nc:
        out = {k: np.array(nc.variables[k][:]) for k in IDENT_KEYS}
        out['time'] = float(nc.variables['matlab_datenum'][...])
        out['UTMzone'] = int(getattr(nc, 'UTMzone', 0))
        slat = np.array(nc.variables['stream_lat'][:])
        slon = np.array(nc.variables['stream_lon'][:])
        starts = np.array(nc.variables['stream_start'][:])
        lens = np.array(nc.variables['stream_len'][:])
        out['stream_eddy'] = np.array(nc.variables['stream_eddy'][:])
        out['stream_lat'] = [slat[s:s + l] for s, l in zip(starts, lens)]
        out['stream_lon'] = [slon[s:s + l] for s, l in zip(starts, lens)]
    return out


TRACK_KEYS = ['lat_center', 'lon_center', 'ellipse_theta', 'eig1', 'eig2',
              'omega', 'direction', 'num_streams']

_TRACK_META = {
    'lat_center': ('degrees_north', ''),
    'lon_center': ('degrees_east', ''),
    'ellipse_theta': ('degrees', 'ellipse rotation'),
    'eig1': ('km', 'ellipse semi-axis 1'),
    'eig2': ('km', 'ellipse semi-axis 2'),
    'omega': ('degree s-1', 'angular velocity'),
    'direction': ('1', '+1 CCW, -1 CW'),
    'num_streams': ('1', 'streamlines (intensity)'),
}


def save_tracks_nc(path, cells, Time, timegap, dist_param, time_param):
    """cells: dict TRACK_KEYS -> list (per track) of per-timestep lists."""
    if not path.endswith('.nc'):
        path = path + '.nc'
    ntracks = len(timegap)
    lens = np.array([len(c) for c in cells['lat_center']], dtype=np.int32)
    npts = int(lens.sum())
    with Dataset(path, 'w', format='NETCDF4') as nc:
        nc.title = 'eddy tracks'
        nc.eddy_track_dist_param = float(dist_param)
        nc.eddy_track_time_param = float(time_param)
        nc.createDimension('track', ntracks)
        nc.createDimension('track_pt', npts)
        vs = nc.createVariable('track_start', 'i4', ('track',))
        vs.long_name = '0-based start offset of each track in the track_pt arrays'
        vl = nc.createVariable('track_len', 'i4', ('track',))
        vl.long_name = 'number of timesteps in each track'
        if ntracks:
            vs[:] = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int32)
            vl[:] = lens
        vg = nc.createVariable('timegap', 'f8', ('track',))
        vg[:] = np.asarray(timegap, dtype=float)

        for k in TRACK_KEYS:
            units, ln = _TRACK_META[k]
            v = nc.createVariable(k, 'f8', ('track_pt',), zlib=True)
            v.units = units
            if ln:
                v.long_name = ln
            if npts:
                v[:] = np.concatenate([np.asarray(c, dtype=float)
                                       for c in cells[k]])
        vt = nc.createVariable('time', 'f8', ('track_pt',), zlib=True)
        vt.units = 'days since 1970-01-01 00:00:00'
        vdn = nc.createVariable('matlab_datenum', 'f8', ('track_pt',), zlib=True)
        if npts:
            tt = np.concatenate([np.asarray(c, dtype=float) for c in Time])
            vt[:] = tt - DATENUM_1970
            vdn[:] = tt


def load_tracks_nc(path):
    """Read tracks -> dict with per-track lists of arrays + timegap/params."""
    with Dataset(path) as nc:
        starts = np.array(nc.variables['track_start'][:])
        lens = np.array(nc.variables['track_len'][:])
        out = {}
        for k in TRACK_KEYS:
            flat = np.array(nc.variables[k][:])
            out[k] = [flat[s:s + l] for s, l in zip(starts, lens)]
        flat = np.array(nc.variables['matlab_datenum'][:])
        out['Time'] = [flat[s:s + l] for s, l in zip(starts, lens)]
        out['timegap'] = np.array(nc.variables['timegap'][:])
        out['eddy_track_dist_param'] = float(nc.eddy_track_dist_param)
        out['eddy_track_time_param'] = float(nc.eddy_track_time_param)
    return out
