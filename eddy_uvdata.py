#!/usr/bin/env python3
"""Python port of eddy_uvdata.m — run the winding-angle eddy identification
on a single netcdf file (or your own lon, lat, u, v data set).

Mirrors the MATLAB example exactly: same bundled Delaware Bay sample, same
subregion, same parameters, same outputs (data/results/*.mat and
data/results_figs/*.png).
"""

import os
import numpy as np
from netCDF4 import Dataset

from eddy_subroutine import eddy_subroutine

os.chdir(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir('data') and os.path.isdir(os.path.join('..', 'data')):
    os.chdir('..')  # repo layout: python/ subfolder, data/ at the repo root

# %% load data
data_dir = 'data/'
fn = '202004030900_hfr_usegc_2km_rtv_uwls_NDBC'  # filename
nc = Dataset(data_dir + fn + '.nc')
lon = np.array(nc.variables['lon'][:], dtype=float)
lat = np.array(nc.variables['lat'][:], dtype=float)
# transpose to (nlon, nlat) to match MATLAB's ncread orientation
u = np.ma.filled(np.squeeze(nc.variables['u'][:]), np.nan).T
v = np.ma.filled(np.squeeze(nc.variables['v'][:]), np.nan).T
time = np.array(nc.variables['time'][:], dtype=float)
nc.close()
lats, lons = np.meshgrid(lat, lon)  # (nlon, nlat), as in MATLAB

params = {}

# %% lon lat limits
params['lonlatlimits'] = 1  # (=1) limit to subregion, 0 for the entire region
params['lat_min'] = 38.5
params['lat_max'] = 39.4
params['lon_min'] = -75.5
params['lon_max'] = -74.5

# %% parameters
params['save_data_name'] = data_dir + 'results/' + fn
params['save_fig_name'] = data_dir + 'results_figs/' + fn

params['plt_debug'] = 1          # debug plots
params['plt_debug_streams'] = 0  # streamline plots
params['save_streams'] = 1       # save lat lon of streamlines
params['plt_quiver'] = 1         # plot quiver
params['plt_debug_sv'] = 1       # save debug plot
params['highrescoast'] = 1       # high resolution coastline
params['bathy'] = 1              # plot bathy levels
params['levels'] = [-10, -25, -50, -100, -250, -500, -1000, -1500,
                    -2000, -2500, -3000, -4000, -5000]  # isobaths
params['plt_close'] = 1          # close figures after each timestep
params['print_res'] = 300        # resolution of plot
params['plt_each_eddy'] = 0      # zooms and plots figure for each eddy
params['plt_min_streams'] = 5

params['winding_thres'] = 300    # how much winding (in degrees) to be an eddy
params['d_thres'] = 10           # streamline start/end max distance in km
params['baddir_thres'] = 15      # breaks at 15 degrees of bad dirs
radar_km_resolution = 2
params['param_center_dist'] = 2 * np.sqrt(2) * radar_km_resolution
params['new_dist_thres'] = 1     # Cahl et al., 2023 clustering method
params['radar_km_resolution'] = radar_km_resolution

params['dij'] = 50  # streamlines in areas of lons[i-dij:i+dij, j-dij:j+dij]

# %% call eddy finding
if __name__ == '__main__':
    eddy_subroutine(params, u, v, time, lons, lats, lon, lat)
