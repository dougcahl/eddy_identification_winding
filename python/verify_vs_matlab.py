#!/usr/bin/env python3
"""Compare the Python pipeline's NetCDF results against MATLAB's .mat results
for every timestep and for the tracks file.

Usage: python verify_vs_matlab.py [path\\to\\matlab\\repo\\data\\results]
"""

import os
import sys
import glob

import numpy as np
from scipy.io import loadmat

from results_io import load_ident_nc, load_tracks_nc

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')

ml_dir = sys.argv[1] if len(sys.argv) > 1 else \
    r'D:\sync\hf\eddy_identification_winding\data\results'

KEYS = ['eddy_center_lat', 'eddy_center_lon', 'eddy_dir', 'eddy_streamlines',
        'eddy_length_x', 'eddy_length_y', 'eddy_ellipse_theta']

# ---- per-timestep identification results ----
py_files = sorted(f for f in glob.glob('data/results/data2_*.nc')
                  if 'tracks' not in f)
n_cmp = 0
count_mismatch = []
max_diff = {k: 0.0 for k in KEYS}
tot_py = 0
tot_ml = 0
for f in py_files:
    base = os.path.splitext(os.path.basename(f))[0]
    mlf = os.path.join(ml_dir, base + '.mat')
    if not os.path.isfile(mlf):
        continue
    a = load_ident_nc(f)
    b = loadmat(mlf)
    n_cmp += 1
    na = len(a['eddy_streamlines'])
    nb = int(np.ravel(b['eddy_streamlines']).size) if b['eddy_streamlines'].size else 0
    tot_py += na
    tot_ml += nb
    if na != nb:
        count_mismatch.append((base, na, nb))
        continue
    for k in KEYS:
        if na:
            d = np.max(np.abs(np.ravel(a[k]) - np.ravel(b[k]).astype(float)))
            max_diff[k] = max(max_diff[k], float(d))

print('timesteps compared: %d' % n_cmp)
print('total eddies: python %d | matlab %d' % (tot_py, tot_ml))
print('eddy-count mismatches: %d %s' % (len(count_mismatch), count_mismatch[:5]))
for k in KEYS:
    print('  max |python - matlab| %-20s %.3e' % (k, max_diff[k]))

# ---- tracks ----
pt = load_tracks_nc('data/results/data2_tracks.nc')
mt = loadmat(os.path.join(ml_dir, 'data2_tracks.mat'))
la = sorted([len(c) for c in pt['lat_center']], reverse=True)
lb = sorted([len(np.ravel(c)) for c in mt['lat_center'][0]], reverse=True)
print('tracks: python %d | matlab %d | lengths equal: %s'
      % (len(pt['timegap']), mt['lat_center'][0].size, la == lb))
dmax = 0.0
for ca, cb in zip(pt['lat_center'], mt['lat_center'][0]):
    ca = np.ravel(ca); cb = np.ravel(cb)
    if len(ca) == len(cb):
        dmax = max(dmax, float(np.max(np.abs(ca - cb))) if len(ca) else 0.0)
    else:
        dmax = np.inf
print('max |track center lat diff|: %.3e deg' % dmax)
print('VERIFY_DONE')
