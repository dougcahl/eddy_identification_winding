#!/usr/bin/env python3
"""Convert the MATLAB .mat data files used by eddy_identification_winding to
NetCDF4 (.nc) — the standard, Python-friendly format for this toolbox.

Usage:
  python convert_data.py data            # data/data2.mat -> data/data2.nc
  python convert_data.py results         # data/results/*.mat -> .nc (ident + tracks)
  python convert_data.py all             # both

Time is stored CF-style as days since 1970-01-01 (time = matlab_datenum -
719529); the original MATLAB datenum is kept in `matlab_datenum` for exact
round-tripping. Variable-length streamline/track series are stored as
concatenated arrays with start/length index vectors (ragged/CSR layout).
"""

import os
import sys
import glob

import numpy as np
from netCDF4 import Dataset
from scipy.io import loadmat

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root

DATENUM_1970 = 719529.0  # MATLAB datenum of 1970-01-01


def convert_input(matfile='data/data2.mat', ncfile='data/data2.nc'):
    m = loadmat(matfile)
    lons = m['lons']; lats = m['lats']
    u1 = m['u1']; v1 = m['v1']
    time1 = m['time1'].ravel()

    with Dataset(ncfile, 'w', format='NETCDF4') as nc:
        nc.title = 'HF radar surface currents (converted from data2.mat)'
        nc.source = 'eddy_identification_winding example data, WERA HF radar'
        nc.history = 'converted from %s' % matfile

        nc.createDimension('i', lons.shape[0])
        nc.createDimension('j', lons.shape[1])
        nc.createDimension('time', len(time1))

        vlon = nc.createVariable('lon', 'f8', ('i', 'j'), zlib=True)
        vlon.units = 'degrees_east'
        vlat = nc.createVariable('lat', 'f8', ('i', 'j'), zlib=True)
        vlat.units = 'degrees_north'
        vt = nc.createVariable('time', 'f8', ('time',))
        vt.units = 'days since 1970-01-01 00:00:00'
        vt.calendar = 'standard'
        vdn = nc.createVariable('matlab_datenum', 'f8', ('time',))
        vdn.long_name = 'original MATLAB datenum (days since year 0)'
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
        vt[:] = time1 - DATENUM_1970
        vdn[:] = time1
        vu[:] = np.moveaxis(u1, 2, 0)
        vv[:] = np.moveaxis(v1, 2, 0)
    print('wrote %s  (%d x %d grid, %d timesteps)'
          % (ncfile, lons.shape[0], lons.shape[1], len(time1)))


def convert_ident_result(matfile, ncfile):
    m = loadmat(matfile)
    n = int(np.ravel(m['eddy_streamlines']).size) if m['eddy_streamlines'].size else 0

    with Dataset(ncfile, 'w', format='NETCDF4') as nc:
        nc.title = 'winding-angle eddy identification results'
        nc.history = 'converted from %s' % os.path.basename(matfile)
        nc.UTMzone = int(np.ravel(m['UTMzone'])[0]) if 'UTMzone' in m else 0

        nc.createDimension('eddy', n)
        vt = nc.createVariable('time', 'f8', ())
        vt.units = 'days since 1970-01-01 00:00:00'
        tval = float(np.ravel(m['time'])[0])
        vt[...] = tval - DATENUM_1970
        vdn = nc.createVariable('matlab_datenum', 'f8', ())
        vdn[...] = tval

        def scal(name, units='', long_name=''):
            v = nc.createVariable(name, 'f8', ('eddy',))
            if units:
                v.units = units
            if long_name:
                v.long_name = long_name
            v[:] = np.ravel(m[name]).astype(float) if n else []
            return v

        scal('eddy_center_lat', 'degrees_north')
        scal('eddy_center_lon', 'degrees_east')
        scal('eddy_dir', '1', '+1 CCW/cyclonic, -1 CW/anticyclonic')
        scal('eddy_angular_vel', 'degree s-1')
        scal('eddy_streamlines', '1', 'number of streamlines in the eddy cluster')
        scal('eddy_length_x', 'km', 'ellipse semi-axis 1')
        scal('eddy_length_y', 'km', 'ellipse semi-axis 2')
        scal('eddy_ellipse_theta', 'degrees', 'ellipse rotation, math notation')

        # ragged streamline storage
        slat = []
        slon = []
        s_eddy = []
        if 'eddy_streamlines_lat' in m and m['eddy_streamlines_lat'].size:
            cl = m['eddy_streamlines_lat']
            cn = m['eddy_streamlines_lon']
            for a in range(cl.shape[0]):
                for b in range(cl.shape[1]):
                    arr = np.ravel(cl[a, b])
                    if arr.size:
                        slat.append(arr)
                        slon.append(np.ravel(cn[a, b]))
                        s_eddy.append(a)
        nstreams = len(slat)
        npts = int(sum(len(s) for s in slat))
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
            v1[:] = np.concatenate(slat)
            v2[:] = np.concatenate(slon)
            v3[:] = np.asarray(s_eddy, dtype=np.int32)
            lens = np.array([len(s) for s in slat], dtype=np.int32)
            v4[:] = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int32)
            v5[:] = lens


def convert_tracks(matfile, ncfile):
    m = loadmat(matfile)
    cells = {k: [np.ravel(c) for c in m[k][0]] for k in
             ['lat_center', 'lon_center', 'ellipse_theta', 'eig1', 'eig2',
              'omega', 'direction', 'num_streams', 'Time']}
    ntracks = len(cells['lat_center'])
    lens = np.array([len(c) for c in cells['lat_center']], dtype=np.int32)
    npts = int(lens.sum())

    with Dataset(ncfile, 'w', format='NETCDF4') as nc:
        nc.title = 'eddy tracks'
        nc.history = 'converted from %s' % os.path.basename(matfile)
        nc.eddy_track_dist_param = float(np.ravel(m['eddy_track_dist_param'])[0])
        nc.eddy_track_time_param = float(np.ravel(m['eddy_track_time_param'])[0])

        nc.createDimension('track', ntracks)
        nc.createDimension('track_pt', npts)
        vs = nc.createVariable('track_start', 'i4', ('track',))
        vs.long_name = '0-based start offset of each track in the track_pt arrays'
        vl = nc.createVariable('track_len', 'i4', ('track',))
        vl.long_name = 'number of timesteps in each track'
        vs[:] = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int32)
        vl[:] = lens
        vg = nc.createVariable('timegap', 'f8', ('track',))
        vg[:] = np.ravel(m['timegap']).astype(float)

        names = {'lat_center': ('degrees_north', ''),
                 'lon_center': ('degrees_east', ''),
                 'ellipse_theta': ('degrees', 'ellipse rotation'),
                 'eig1': ('km', 'ellipse semi-axis 1'),
                 'eig2': ('km', 'ellipse semi-axis 2'),
                 'omega': ('degree s-1', 'angular velocity'),
                 'direction': ('1', '+1 CCW, -1 CW'),
                 'num_streams': ('1', 'streamlines (intensity)')}
        for k, (units, ln) in names.items():
            v = nc.createVariable(k, 'f8', ('track_pt',), zlib=True)
            v.units = units
            if ln:
                v.long_name = ln
            if npts:
                v[:] = np.concatenate(cells[k]).astype(float)
        vt = nc.createVariable('time', 'f8', ('track_pt',), zlib=True)
        vt.units = 'days since 1970-01-01 00:00:00'
        vdn = nc.createVariable('matlab_datenum', 'f8', ('track_pt',), zlib=True)
        if npts:
            tt = np.concatenate(cells['Time']).astype(float)
            vt[:] = tt - DATENUM_1970
            vdn[:] = tt
    print('wrote %s  (%d tracks, %d track points)' % (ncfile, ntracks, npts))


def convert_results():
    files = sorted(glob.glob('data/results/*.mat'))
    count = 0
    for f in files:
        out = f[:-4] + '.nc'
        if 'tracks' in os.path.basename(f):
            convert_tracks(f, out)
        else:
            convert_ident_result(f, out)
            count += 1
    print('converted %d identification result files to .nc' % count)


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('data', 'all'):
        convert_input()
    if what in ('results', 'all'):
        convert_results()
