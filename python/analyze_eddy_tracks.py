#!/usr/bin/env python3
"""Python port of analyze_eddy_tracks.m — plot the eddy tracks found by
eddy_tracking.py and compute track statistics.

Figure 1: track map (blue = CCW/cyclonic, red = CW/anticyclonic; green o =
start, black x = end) over coastline + bathymetry.
Figure 2: histograms of duration, intensity, speed, and diameters.
Both are saved to data/results_figs/.
"""

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

from utm_dc import geog2utm_nodisp
from etopo1_reader import read_etopo1

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root

plt_individual = 0
highrescoast = 1
lw = 1
min_steps = 12   # minimum track length to plot
radar_dt = 0.5   # timestep in hours

# %% load
name_dir = 'data/results/'
name_pre = 'data2_'
t = loadmat(name_dir + name_pre + 'tracks.mat')
data = loadmat('data/data2.mat')
lons = data['lons']; lats = data['lats']

lon_center = [np.ravel(c) for c in t['lon_center'][0]]
lat_center = [np.ravel(c) for c in t['lat_center'][0]]
omega = [np.ravel(c) for c in t['omega'][0]]
num_streams = [np.ravel(c) for c in t['num_streams'][0]]
eig1 = [np.ravel(c) for c in t['eig1'][0]]
eig2 = [np.ravel(c) for c in t['eig2'][0]]
direction = [np.ravel(c) for c in t['direction'][0]]
Time = [np.ravel(c) for c in t['Time'][0]]
timegap = np.ravel(t['timegap'])
eddy_track_time_param = float(np.ravel(t['eddy_track_time_param'])[0])

# %% setup
lon_min = float(np.floor(10 * np.nanmin(lons)) / 10)
lon_max = float(np.ceil(10 * np.nanmax(lons)) / 10)
lat_min = float(np.floor(10 * np.nanmin(lats)) / 10)
lat_max = float(np.ceil(10 * np.nanmax(lats)) / 10)

# %% per-track statistics
ntracks = len(timegap)
print('max gap between eddy identification in timesteps = %g' % eddy_track_time_param)
print('number of tracks = %d' % ntracks)

eddy_track_length = np.array([len(d) for d in direction])
intensity = np.array([np.mean(s) for s in num_streams])
diameter = np.array([np.mean(e1 + e2) for e1, e2 in zip(eig1, eig2)])
maxdiameter = np.array([np.max(e1 + e2) for e1, e2 in zip(eig1, eig2)])
mindiameter = np.array([np.min(e1 + e2) for e1, e2 in zip(eig1, eig2)])
dist = np.full(ntracks, np.nan)
speed = np.full(ntracks, np.nan)
for i in range(ntracks):
    if eddy_track_length[i] > 1:
        x, y, _ = geog2utm_nodisp(lon_center[i], lat_center[i])
        d1 = np.hypot(x[-1] - x[0], y[-1] - y[0]) * 1000
        dt = (Time[i][-1] - Time[i][0]) * 24 * 3600
        dist[i] = d1
        speed[i] = d1 / dt

print('max eddy track length in timesteps = %d' % np.max(eddy_track_length))

# %% track map
order = np.argsort(-eddy_track_length)  # descending

fig1 = plt.figure(figsize=(11, 9), facecolor='w')
ax = fig1.add_subplot(111)
try:
    elev, elon, elat = read_etopo1(lon_min - 0.05, lon_max + 0.05,
                                   lat_min - 0.05, lat_max + 0.05)
    levels = [-2500, -1000, -500, -250, -100, -50, -25, -10]
    ax.contour(elon, elat, elev, levels=levels, colors=[(0.5, 0.5, 0.5)],
               linewidths=0.5, linestyles='-', zorder=1)
    land = np.ma.masked_less(elev, 0)
    ax.contourf(elon, elat, land, levels=[0, elev.max() + 1],
                colors=[(0.8, 0.8, 0.8)], zorder=2)
    ax.contour(elon, elat, elev, levels=[0], colors='k', linewidths=0.7, zorder=3)
except FileNotFoundError as e:
    print('warning: %s (skipping bathymetry/coastline)' % e)

counter = 0
title_set = False
for i in order:
    if eddy_track_length[i] < min_steps:
        continue
    counter += 1
    w1 = np.nanmean(omega[i])
    w2 = w1 * 3600 * 24 / 360  # rotations per day
    lonp = lon_center[i]
    latp = lat_center[i]
    color = 'b' if w1 > 0 else 'r'
    ax.plot(lonp, latp, '-' + color + '.', linewidth=lw, markersize=5, zorder=4)
    ax.plot(lonp[0], latp[0], 'o', linewidth=2, markerfacecolor='none',
            markeredgecolor=(0, 0.8, 0), markeredgewidth=2, zorder=5)
    ax.plot(lonp[-1], latp[-1], 'kx', markersize=9, markeredgewidth=2, zorder=5)
    if not title_set:
        ax.set_title('longest track = %d timesteps  ;   %.1f rot/day'
                     % (eddy_track_length[i], round(w2, 1)))
        title_set = True

ax.set_xlim(lon_min, lon_max)
ax.set_ylim(lat_min, lat_max)
ax.set_aspect(1 / np.cos(np.deg2rad((lat_min + lat_max) / 2)))
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
fig1.savefig('data/results_figs/data2_analysis_tracks.png', dpi=200,
             bbox_inches='tight')
plt.close(fig1)
print('%d tracks with >= %d timesteps plotted' % (counter, min_steps))

# %% histogram plots
fig2 = plt.figure(figsize=(13, 7), facecolor='w')
panels = [
    (eddy_track_length * radar_dt, 'Eddy duration (hours)'),
    (intensity, 'intensity (avg streamlines)'),
    (speed[~np.isnan(speed)], 'speed (m/s)'),
    (diameter, 'avg diameter (km)'),
    (mindiameter, 'min diameter (km)'),
    (maxdiameter, 'max diameter (km)'),
]
for ip, (vals, lab) in enumerate(panels, start=1):
    axh = fig2.add_subplot(2, 3, ip)
    axh.hist(vals, bins=10, edgecolor='k', linewidth=0.5)
    axh.set_xlabel(lab)
    if ip == 1:
        axh.set_ylabel('Number of eddies')
    axh.grid(True, alpha=0.4)
fig2.tight_layout()
fig2.savefig('data/results_figs/data2_analysis_hist.png', dpi=200,
             bbox_inches='tight')
plt.close(fig2)
print('figures saved to data/results_figs/')
