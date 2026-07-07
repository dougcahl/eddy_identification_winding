#!/usr/bin/env python3
"""Python port of eddy_tracking.m — link per-timestep eddy identifications
into tracks.

An eddy at this timestep continues a track when (Cahl et al., 2023) the new
center falls inside the previous eddy's fitted ellipse (inpolygon test), or
the centers are within eddy_track_dist_param km; a track can bridge up to
eddy_track_time_param timesteps without an identification.

Saves data/results/data2_tracks.mat with the same variables as MATLAB.
"""

import os
import glob

import numpy as np
from matplotlib.path import Path

from utm_dc import geog2utm_nodisp
from results_io import load_ident_nc, save_tracks_nc

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root

# %% parameters
name_dir = 'data/results/'
name_pre = 'data2_'
outfile = name_dir + name_pre + 'tracks'

eddy_track_dist_param = 10  # km from previous eddy center to be a continuation
eddy_track_time_param = 2   # max timesteps between identifications
new_dist_thres = 1


def _inpoly(px, py, polyx, polyy):
    polyx = np.atleast_1d(polyx)
    polyy = np.atleast_1d(polyy)
    if polyx.size < 3:  # degenerate polygon: MATLAB inpolygon returns false
        return False    # (the in2 test below builds a 1-point 'ellipse', as in MATLAB)
    p = Path(np.column_stack([polyx, polyy]))
    return bool(p.contains_point((px, py)))


# %% load data
files = []
for fn in sorted(glob.glob(name_dir + name_pre + '*.nc')):
    base = os.path.basename(fn)
    if 'tracks' in base:
        continue
    files.append(fn)

recs = [load_ident_nc(fn) for fn in files]

# %% sort by time
recs.sort(key=lambda r: r['time'])
time = np.array([r['time'] for r in recs])
eddy_center_lat = [r['eddy_center_lat'] for r in recs]
eddy_center_lon = [r['eddy_center_lon'] for r in recs]
eddy_dir = [r['eddy_dir'] for r in recs]
eddy_ellipse_theta = [r['eddy_ellipse_theta'] for r in recs]
eddy_length_x = [r['eddy_length_x'] for r in recs]
eddy_length_y = [r['eddy_length_y'] for r in recs]
eddy_angular_vel = [r['eddy_angular_vel'] for r in recs]
eddy_streamlines = [r['eddy_streamlines'] for r in recs]

total_eddies = sum(len(d) for d in eddy_dir)
print('Total eddies = %d' % total_eddies)
print('Max eddies in a single timestep = %d' % max(len(d) for d in eddy_dir))
print('Max streamlines for an eddy = %d'
      % max((int(np.max(s)) if s.size else 0) for s in eddy_streamlines))

# %% no tracking for first timestep
it = 0
i_cluster = len(eddy_center_lat[it])
lat_center = []; lon_center = []; ellipse_theta = []
eig1 = []; eig2 = []; omega = []; direction = []; num_streams = []
Time = []; timegap = []
for i in range(i_cluster):
    lat_center.append([eddy_center_lat[it][i]])
    lon_center.append([eddy_center_lon[it][i]])
    ellipse_theta.append([eddy_ellipse_theta[it][i]])
    eig1.append([eddy_length_x[it][i]])
    eig2.append([eddy_length_y[it][i]])
    omega.append([eddy_angular_vel[it][i]])
    direction.append([eddy_dir[it][i]])
    num_streams.append([eddy_streamlines[it][i]])
    Time.append([time[it]])
    timegap.append(1)

total_tracks = i_cluster
print('Starting track finding routine')
print('Total number of timestep to analyze = %d' % len(time))

i_ellipse = np.arange(30, 361, 30)
cosI = np.cos(np.deg2rad(i_ellipse))
sinI = np.sin(np.deg2rad(i_ellipse))

for it in range(1, len(time)):
    i_cluster = len(eddy_center_lat[it])
    i1_start = 0
    for i in range(i_cluster):
        curr_lat_center = eddy_center_lat[it][i]
        if np.isnan(curr_lat_center):
            continue
        curr_lon_center = eddy_center_lon[it][i]
        curr_ellipse_x = eddy_length_x[it][i]
        curr_ellipse_y = eddy_length_y[it][i]
        curr_ellipse_th = eddy_ellipse_theta[it][i]

        succ = 0
        alive = [k for k in range(total_tracks) if timegap[k] <= eddy_track_time_param]
        i1_start = alive[0] if alive else total_tracks
        for i1 in range(i1_start, total_tracks):
            eddy_continuation = 0
            if new_dist_thres:  # Cahl et al., 2023
                # previous eddy's ellipse (in km), centered on the offset of the
                # old center relative to the current center
                ellipse_x = eig1[i1][-1] * cosI
                ellipse_y = eig2[i1][-1] * sinI
                th = np.deg2rad(ellipse_theta[i1][-1])
                ex = ellipse_x * np.cos(th) - ellipse_y * np.sin(th)
                ey = ellipse_x * np.sin(th) + ellipse_y * np.cos(th)

                x, y, _ = geog2utm_nodisp(lon_center[i1][-1], lat_center[i1][-1],
                                          curr_lon_center, curr_lat_center)
                x = float(x[0]); y = float(y[0])
                ex = ex + x
                ey = ey + y

                # is the current center within the old eddy's ellipse
                in1 = _inpoly(0.0, 0.0, ex, ey)

                # (as in MATLAB: this builds a single point, so in2 is
                # effectively always false — kept for faithfulness)
                th = np.deg2rad(curr_ellipse_th)
                ex1 = curr_ellipse_x * np.cos(th) - curr_ellipse_y * np.sin(th)
                ey1 = curr_ellipse_x * np.sin(th) + curr_ellipse_y * np.cos(th)
                in2 = _inpoly(x, y, ex1, ey1)

                if in1 or in2:
                    eddy_continuation = 1
                else:
                    if np.hypot(x, y) < eddy_track_dist_param:
                        eddy_continuation = 1
            else:  # Sadarjoen and Post 2000
                x, y, _ = geog2utm_nodisp(lon_center[i1][-1], lat_center[i1][-1],
                                          curr_lon_center, curr_lat_center)
                if np.hypot(float(x[0]), float(y[0])) < eddy_track_dist_param:
                    eddy_continuation = 1

            if eddy_continuation:
                for i_tgap in range(1, eddy_track_time_param + 1):
                    if timegap[i1] <= i_tgap:
                        lat_center[i1].append(eddy_center_lat[it][i])
                        lon_center[i1].append(eddy_center_lon[it][i])
                        ellipse_theta[i1].append(eddy_ellipse_theta[it][i])
                        eig1[i1].append(eddy_length_x[it][i])
                        eig2[i1].append(eddy_length_y[it][i])
                        omega[i1].append(eddy_angular_vel[it][i])
                        direction[i1].append(eddy_dir[it][i])
                        num_streams[i1].append(eddy_streamlines[it][i])
                        Time[i1].append(time[it])
                        timegap[i1] = 0
                        succ = 1
                        break

        if succ == 0:  # create a new track
            total_tracks += 1
            lat_center.append([eddy_center_lat[it][i]])
            lon_center.append([eddy_center_lon[it][i]])
            ellipse_theta.append([eddy_ellipse_theta[it][i]])
            eig1.append([eddy_length_x[it][i]])
            eig2.append([eddy_length_y[it][i]])
            omega.append([eddy_angular_vel[it][i]])
            direction.append([eddy_dir[it][i]])
            num_streams.append([eddy_streamlines[it][i]])
            Time.append([time[it]])
            timegap.append(0)

    for i1 in range(i1_start, total_tracks):
        timegap[i1] += 1

print('Finished eddy track finding')
print('number of tracks = %d' % total_tracks)

# %% save (native NetCDF)
save_tracks_nc(outfile,
               {'lat_center': lat_center, 'lon_center': lon_center,
                'ellipse_theta': ellipse_theta, 'eig1': eig1, 'eig2': eig2,
                'omega': omega, 'direction': direction,
                'num_streams': num_streams},
               Time, timegap, eddy_track_dist_param, eddy_track_time_param)
print('saved %s.nc' % outfile)
