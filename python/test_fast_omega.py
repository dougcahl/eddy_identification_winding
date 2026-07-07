#!/usr/bin/env python3
"""A/B test: run eddy_subroutine with fast_omega=1 vs fast_omega=0 on several
timesteps and report exactly which saved quantities differ."""

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root

STEPS = [1, 25, 50, 61, 75, 101]  # 1-based


def run_pair(i):
    from scipy.io import loadmat
    from eddy_uvdata_loop import load_data, make_params
    from eddy_subroutine import eddy_subroutine
    lons, lats, u1, v1, time1 = load_data()
    out = {}
    for fast in (1, 0):
        params = make_params()
        params['fast_omega'] = fast
        params['save_data_name'] = 'data/testing/ab_%d_fast%d' % (i, fast)
        params['save_fig_name'] = 'data/testing/ab_%d_fast%d' % (i, fast)
        eddy_subroutine(params, u1[:, :, i - 1], v1[:, :, i - 1], time1[i - 1],
                        lons, lats)
        out[fast] = loadmat(params['save_data_name'] + '.mat')
    a, b = out[1], out[0]
    geom_keys = ['eddy_center_lat', 'eddy_center_lon', 'eddy_dir',
                 'eddy_streamlines', 'eddy_length_x', 'eddy_length_y',
                 'eddy_ellipse_theta', 'eddy_streamlines_lat',
                 'eddy_streamlines_lon']

    def eq(x, y):
        if x.dtype == object:
            return x.shape == y.shape and all(
                np.array_equal(np.ravel(x[i2, j2]), np.ravel(y[i2, j2]))
                for i2 in range(x.shape[0]) for j2 in range(x.shape[1]))
        return np.array_equal(x, y)

    geom_same = all(eq(a[k], b[k]) for k in geom_keys)
    oa = np.ravel(a['eddy_angular_vel'])
    ob = np.ravel(b['eddy_angular_vel'])
    nan_a, nan_b = np.isnan(oa), np.isnan(ob)
    both = ~nan_a & ~nan_b
    if both.any():
        rel = float(np.max(np.abs(oa[both] - ob[both]) / np.abs(ob[both])))
    else:
        rel = 0.0
    return (i, len(oa), geom_same, int(nan_a.sum()), int(nan_b.sum()),
            bool(np.array_equal(oa, ob, equal_nan=True)), rel)


if __name__ == '__main__':
    os.makedirs('data/testing', exist_ok=True)
    with ProcessPoolExecutor(max_workers=len(STEPS)) as ex:
        results = list(ex.map(run_pair, STEPS))
    print('%-5s %-7s %-18s %-12s %-12s %-12s %s' %
          ('step', 'eddies', 'geometry_identical', 'omega_nan_f', 'omega_nan_x',
           'omega_equal', 'max_rel_omega_diff'))
    for r in results:
        print('%-5d %-7d %-18s %-12d %-12d %-12s %.2e' %
              (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))
