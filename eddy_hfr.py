#!/usr/bin/env python3
"""Download HF radar surface currents from ERDDAP and run the full
winding-angle eddy pipeline: identification (parallel over CPU cores),
tracking, and track analysis, with a small set of overview figures.

Example:
  python eddy_hfr.py --dataset ucsdHfrW6 --start 2026-07-25 --end 2026-07-27 \
                     --bbox -125 -120.5 33.5 38.5 --min-steps 6

Outputs (prefix = hfr_<dataset>_<YYYYMMDD of start>, or --prefix):
  data/<prefix>.nc                      downloaded input currents
  data/results/<prefix>_<i>.nc          per-timestep identifications
  data/results/<prefix>_tracks.nc       eddy tracks
  data/results_figs/<prefix>_*.png      figures
"""

import argparse
import os
import sys
import time as _time
from datetime import datetime

import numpy as np
from netCDF4 import Dataset

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root

from hfr_download import download_hfr, resolve_dataset
from eddy_tracking import run_tracking
from analyze_eddy_tracks import run_analysis
from results_io import load_ident_nc


def make_hfr_params(res_km, plt_debug=0):
    """The repo's standard identification parameters, resolution-scaled so
    the physical streamline search window stays ~200 km."""
    params = {}
    params['lonlatlimits'] = 0

    params['plt_debug'] = plt_debug
    params['plt_debug_streams'] = 0
    params['save_streams'] = 1
    params['plt_quiver'] = 1
    params['plt_debug_sv'] = 1
    params['highrescoast'] = 1
    params['bathy'] = 1
    params['levels'] = [-10, -25, -50, -100, -250, -500, -1000, -1500,
                        -2000, -2500, -3000, -4000, -5000]
    params['plt_close'] = 1
    params['print_res'] = 150
    params['plt_each_eddy'] = 0
    params['plt_min_streams'] = 5

    params['winding_thres'] = 300
    params['d_thres'] = 10
    params['baddir_thres'] = 15
    params['param_center_dist'] = 2 * np.sqrt(2) * res_km
    params['new_dist_thres'] = 1
    params['radar_km_resolution'] = res_km
    params['dij'] = int(round(100 / res_km))
    return params


def load_hfr_nc(path):
    """Load a toolbox-standard input file written by hfr_download (axes are
    already ascending; do NOT apply eddy_uvdata_loop.load_data, whose flip
    is specific to the data2 example)."""
    with Dataset(path) as nc:
        lons = np.array(nc.variables['lon'][:])
        lats = np.array(nc.variables['lat'][:])
        u1 = np.moveaxis(np.ma.filled(nc.variables['u'][:], np.nan), 0, 2)
        v1 = np.moveaxis(np.ma.filled(nc.variables['v'][:], np.nan), 0, 2)
        time1 = np.array(nc.variables['matlab_datenum'][:])
    assert lats[0, 1] > lats[0, 0], 'latitude axis must ascend'
    assert lons[1, 0] > lons[0, 0], 'longitude axis must ascend'
    return lons, lats, u1, v1, time1


def run_one(i, nc_path, prefix, res_km, plt_debug):
    from eddy_subroutine import eddy_subroutine
    lons, lats, u1, v1, time1 = load_hfr_nc(nc_path)
    params = make_hfr_params(res_km, plt_debug)
    params['save_data_name'] = 'data/results/%s_%d' % (prefix, i + 1)
    params['save_fig_name'] = 'data/results_figs/%s_%d' % (prefix, i + 1)
    t0 = _time.time()
    eddy_subroutine(params, u1[:, :, i], v1[:, :, i], time1[i], lons, lats)
    return i + 1, _time.time() - t0


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--dataset', required=True,
                   help='ERDDAP id (ucsdHfrE1/E2/E6/W1/W2/W6) or preset '
                        '(USEGC/USWC with --res)')
    p.add_argument('--res', type=int, default=None, choices=(1, 2, 6))
    p.add_argument('--start', required=True, help='YYYY-MM-DD[THH:MM]')
    p.add_argument('--end', required=True)
    p.add_argument('--bbox', nargs=4, type=float, required=True,
                   metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'))
    p.add_argument('--max-hdop', type=float, default=None,
                   help='QC: mask vectors with hdop above this (e.g. 1.25)')
    p.add_argument('--workers', type=int, default=None)
    p.add_argument('--figs', choices=('quick', 'all', 'none'), default='quick',
                   help='quick: first/last/busiest detection maps + analysis '
                        'figures (default); all: every timestep; none')
    p.add_argument('--min-steps', type=int, default=6,
                   help='min track length (timesteps) to draw on the map')
    p.add_argument('--track-dist', type=float, default=10,
                   help='tracking continuation distance (km)')
    p.add_argument('--track-gap', type=int, default=2,
                   help='max TIMESTEPS bridged between identifications')
    p.add_argument('--skip-download', action='store_true')
    p.add_argument('--prefix', default=None)
    args = p.parse_args()

    dataset_id, res_km = resolve_dataset(args.dataset, args.res)
    if res_km == 1:
        print('note: 1 km datasets are 10-100x slower to process than 6 km; '
              'keep the bbox small')
    start_tag = args.start[:10].replace('-', '')
    prefix = args.prefix or 'hfr_%s_%s' % (dataset_id, start_tag)
    nc_path = 'data/%s.nc' % prefix

    # ---------------- download ----------------
    if args.skip_download and os.path.isfile(nc_path):
        print('using existing %s' % nc_path)
    else:
        download_hfr(dataset_id, args.start, args.end, args.bbox, nc_path,
                     max_hdop=args.max_hdop)

    lons, lats, u1, v1, time1 = load_hfr_nc(nc_path)
    nt = len(time1)
    dt_hours = float(np.median(np.diff(time1)) * 24) if nt > 1 else 1.0
    print('grid %dx%d, %d timesteps, dt=%.2f h, %.0f%% mean valid cells'
          % (u1.shape[0], u1.shape[1], nt, dt_hours,
             100 * np.isfinite(u1).mean()))

    # ---------------- identification (parallel) ----------------
    from concurrent.futures import ProcessPoolExecutor, as_completed
    plt_all = 1 if args.figs == 'all' else 0
    fig_steps = set()
    if args.figs == 'quick':
        fig_steps = {0, nt - 1}
    workers = args.workers or max(1, (os.cpu_count() or 4) - 2)
    print('identifying eddies in %d timesteps on %d workers' % (nt, workers))
    t0 = _time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_one, i, nc_path, prefix, res_km,
                          plt_all or (1 if i in fig_steps else 0))
                for i in range(nt)]
        done = 0
        for f in as_completed(futs):
            i1, dt = f.result()
            done += 1
            if done % 10 == 0 or done == nt:
                print('%d/%d timesteps done' % (done, nt), flush=True)
    print('identification: %.1f s' % (_time.time() - t0))

    # eddy counts; busiest-timestep figure for --figs quick
    counts = []
    for i in range(nt):
        r = load_ident_nc('data/results/%s_%d.nc' % (prefix, i + 1))
        counts.append(len(r['eddy_streamlines']))
    print('total eddies identified: %d (max %d in one timestep)'
          % (sum(counts), max(counts)))
    if args.figs == 'quick' and counts:
        busiest = int(np.argmax(counts))
        if busiest not in fig_steps:
            run_one(busiest, nc_path, prefix, res_km, 1)
            print('busiest timestep %d re-rendered (%d eddies)'
                  % (busiest + 1, counts[busiest]))

    # ---------------- tracking + analysis ----------------
    run_tracking(name_pre=prefix + '_',
                 eddy_track_dist_param=args.track_dist,
                 eddy_track_time_param=args.track_gap)
    figs = []
    if args.figs != 'none' or True:  # analysis figures are always cheap
        figs = run_analysis(name_pre=prefix + '_', input_nc=nc_path,
                            radar_dt=dt_hours, min_steps=args.min_steps)

    print('\nDONE. prefix=%s' % prefix)
    print('input:   %s' % nc_path)
    print('results: data/results/%s_*.nc (+ %s_tracks.nc)' % (prefix, prefix))
    print('figures: data/results_figs/%s_*' % prefix)
    return 0


if __name__ == '__main__':
    sys.exit(main())
