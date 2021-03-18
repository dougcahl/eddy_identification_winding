#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 5 12:03:02 2021

@author: doug
"""


import datetime
import os
import ftplib
import time
import math
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import utm
import sys  
from streamplot_dc import streamplot_dc # custom streamplot function that only calculates a single stream line
import csv
import urllib.request

sleeptime = 10 # sleep between analysis in minutes

debug_plot = 0 # change to 0 when deployed



while True:
#    for i_km in range(1,3):
#        if i_km == 0:
#            radar_km_resolution = 6         #% 1, 2 or 6 km nc file
#        if i_km == 1:
#            radar_km_resolution = 2         #% 1, 2 or 6 km nc file
#        if i_km == 2:
#            radar_km_resolution = 1         #% 1, 2 or 6 km nc file

  
    avg_25hr                = 1         # set to 0 to use hourly currents. set to 1 for 25 hour average currents
    radar_km_resolution     = 6         # 1, 2 or 6 km nc file, does not support west coast 500m
    eddy_track_dist_param   = 25        # 4.0; % distance in km for eddy to be continuation
    eddy_track_time_param   = 5         # how many hours can there be between identified eddies for it to be considered the same eddy
    debug_info              = 1         # print when eddy is found for each streamline
    res_x_mult              = 1         # for higher resolution, if >1 results in this squared processing time (i.e. 3 => 9 times longer)
    winding_thres           = 270       # how much winding (in degrees) to be a eddy
    d_thres                 = 25        # streamline starting and ending max distance in km
    baddir_thres            = 15        # breaks at 15 degrees of bad dirs
    param_center_dist       = 25        # how many km are clusters together
    min_pts                 = 5         # min number of grid points streamline has
    grid_deg_size           = 2         # utm grid size for streamlines (ex. 1x1 grid is about 100km x 100km )
  
    ### region selection
    region = 'USEGC'    # 1,2,6 km
    fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'

#    region = 'GAK'      # 1,2,6 km
#    fullregionname = 'HFRADAR_US_Gulf_of_Alaska'
#
#    region = 'AKNS'     # 6 km
#    fullregionname = 'HFRADAR_Alaska_-_North_Slope'
#
#    region = 'USWC'     # 500m 1,2,6 km (500m is not supported)
#    fullregionname = 'HFRADAR_US_West_Coast'
#
#    region = 'USHI'     # 1,2,6 km
#    fullregionname = 'HFRADAR_US_Hawaii'
#
#    region = 'PRVI'     # 2,6 km
#    fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands'



    # directory where hf currents are saved
    din = 'Data/analysis_data/' # fiename in hours since 2019-09-29 20:00:00.000 UTC
    if not os.path.isdir(din):
        os.mkdir(din) # create directory
    din = din + 'km' + str(radar_km_resolution) + '/'
    if not os.path.isdir(din):
        os.mkdir(din) # create directory

    # results directory for ftp upload
    fdir = 'Data/eddy_tracks/'
    if not os.path.isdir(fdir):
        os.mkdir(fdir) # create directory
    fdir = fdir + 'km' + str(radar_km_resolution) + '/'
    if not os.path.isdir(fdir):
        os.mkdir(fdir) # create directory

    file1 = 'eddy_tracks_' + str(radar_km_resolution) + 'km.csv'
    file2 = 'eddy_tracks_' + str(radar_km_resolution) + 'km.npz'
    
    now = datetime.datetime.utcnow()
    start = datetime.datetime(2019,9,29)
    diff = now - start 
    days, seconds = diff.days, diff.seconds
    hours_now = days * 24 + seconds // 3600
    hours_now = hours_now - 20; # change to hours_now = hours_now - 20-100;
    

    files = sorted(os.listdir(din))
    if len(files) == 0:
        numfiles = 1
        nc_time_start = hours_now - 20
    else:
        lastfile = files[-1]
        lasthr = int(lastfile[-8:-3]) # most recent file
        nc_time_start = lasthr + 1
        numfiles = hours_now - lasthr
    
    print('current time = ' + str(hours_now))
    print('25 hr average 12.5 hours behind')
    print('last file +1 = ' + str(nc_time_start))
    print('catchup num  = ' + str(numfiles))
#    numfiles = 1


    exec(open('eddy1_identification.py').read()) # does eddy analysis and debug plots
    exec(open('eddy2_track_eddies.py').read()) # tracks the eddys through time
    exec(open('eddy3_analyze_eddy_tracks.py').read()) # output of last timestep to csv for web graphics
#    break
    
    # uploads results to server
    print('sending results')
    # try:
        # session = ftplib.FTP('i','ftp_user','ftp_pass') # ftp server
        # file = open(fdir + file1,'rb')                  # file to send
        # session.storbinary('STOR ' + file1, file)       # send the file
        # file = open(fdir + file2,'rb')                 
        # session.storbinary('STOR ' + file2, file)     
        # file.close()                                    # close file and FTPsession = ftplib.FTP('104.128.235.148','ftp_home','homehome') 
        # session.quit()
    # except:
        # print('ftp server down, no results going to webserver')
    
        
    print('sleep for ' + str(sleeptime) + ' minutes ...') 
    time.sleep(60*sleeptime) # sleep until next run



