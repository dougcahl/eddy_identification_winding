#!/usr/bin/env python3
"""Python port of eddy_uvdata_loop.m — run eddy identification on every
timestep of a lon, lat, u, v, time data set (data/data2.mat), in parallel
across CPU cores (each timestep is independent).

Outputs one data/results/data2_<i>.mat per timestep, same as MATLAB.
"""

import os
import sys
import time as _time

import numpy as np
from netCDF4 import Dataset

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root


def make_params():
    params = {}
    params['lonlatlimits'] = 0  # 0 = entire region

    params['plt_debug'] = 0          # debug plots
    params['plt_debug_streams'] = 0
    params['save_streams'] = 1
    params['plt_quiver'] = 1
    params['plt_debug_sv'] = 1
    params['highrescoast'] = 1
    params['bathy'] = 1
    params['levels'] = [-10, -25, -50, -100, -250, -500, -1000, -1500,
                        -2000, -2500, -3000, -4000, -5000]
    params['plt_close'] = 1
    params['print_res'] = 300
    params['plt_each_eddy'] = 0
    params['plt_min_streams'] = 5

    params['winding_thres'] = 300
    params['d_thres'] = 10
    params['baddir_thres'] = 15
    radar_km_resolution = 3
    params['param_center_dist'] = 2 * np.sqrt(2) * radar_km_resolution
    params['new_dist_thres'] = 1
    params['radar_km_resolution'] = radar_km_resolution
    params['dij'] = 50
    return params


def load_data():
    # data2.nc is the converted data2.mat (see convert_data.py)
    with Dataset('data/data2.nc') as nc:
        lats = np.array(nc.variables['lat'][:])
        lons = np.array(nc.variables['lon'][:])
        u1 = np.moveaxis(np.ma.filled(nc.variables['u'][:], np.nan), 0, 2)
        v1 = np.moveaxis(np.ma.filled(nc.variables['v'][:], np.nan), 0, 2)
        # MATLAB datenum, so saved results stay comparable with the MATLAB repo
        time1 = np.array(nc.variables['matlab_datenum'][:])
    # WERA site: flip so lats increase with increasing index (as in MATLAB)
    lats = np.fliplr(lats)
    lons = np.fliplr(lons)
    u1 = np.flip(u1, axis=1)
    v1 = np.flip(v1, axis=1)
    return lons, lats, u1, v1, time1


def run_one(i):
    """Run identification for timestep i (0-based); files named 1-based."""
    from eddy_subroutine import eddy_subroutine
    lons, lats, u1, v1, time1 = load_data()
    params = make_params()
    params['save_data_name'] = 'data/results/data2_%d' % (i + 1)
    params['save_fig_name'] = 'data/results_figs/data2_%d' % (i + 1)
    t0 = _time.time()
    eddy_subroutine(params, u1[:, :, i], v1[:, :, i], time1[i], lons, lats)
    return i + 1, _time.time() - t0


if __name__ == '__main__':
    from concurrent.futures import ProcessPoolExecutor, as_completed

    lons, lats, u1, v1, time1 = load_data()
    n = len(time1)
    idx = range(n)
    if len(sys.argv) > 1:  # optional: run a single timestep, e.g. python eddy_uvdata_loop.py 61
        idx = [int(sys.argv[1]) - 1]

    workers = max(1, (os.cpu_count() or 4) - 2)
    print('running %d timesteps on %d workers' % (len(list(idx)), workers))
    t0 = _time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_one, i) for i in idx]
        done = 0
        for f in as_completed(futs):
            i1, dt = f.result()
            done += 1
            print('timestep %d done in %.1f s  (%d/%d)' % (i1, dt, done, len(futs)),
                  flush=True)
    print('LOOP_DONE %.1f s' % (_time.time() - t0))
