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

plt.rcParams.update({'figure.max_open_warning': 50})
sleeptime = 10 # sleep between analysis in minutes

debug_plot = 1 # change to 0 when deployed



while True:
    for i_km in range(0,3):
        if i_km == 0:
            radar_km_resolution = 6         #% 1, 2 or 6 km nc file
        if i_km == 1:
            radar_km_resolution = 2         #% 1, 2 or 6 km nc file
        if i_km == 2:
            radar_km_resolution = 1         #% 1, 2 or 6 km nc file

    
    #    radar_km_resolution = 6         #% 1, 2 or 6 km nc file


        # directory where hf currents are saved
        din = 'Data/analysis_data/' # fiename in hours since 2019-09-29 20:00:00.000 UTC
        din = din + 'km' + str(radar_km_resolution) + '/'
        if not os.path.isdir(din):
            os.mkdir(din) # create directory

        # results directory for ftp upload
        fdir = 'Data/eddy_tracks/'
        fdir = fdir + 'km' + str(radar_km_resolution) + '/'
        if not os.path.isdir(fdir):
            os.mkdir(fdir) # create directory

        file1 = 'eddy_tracks_' + str(radar_km_resolution) + 'km.csv'
        file2 = 'eddy_tracks_' + str(radar_km_resolution) + 'km.npz'
        
        now = datetime.datetime.now()
        start = datetime.datetime(2019,9,29)
        diff = now - start 
        days, seconds = diff.days, diff.seconds
        hours_now = days * 24 + seconds // 3600
        hours_now = hours_now - 20; # change to hours_now = hours_now - 20-100;
        print(hours_now)

        files = sorted(os.listdir(din))
        if len(files) == 0:
            numfiles = 1
            nc_time_start = hours_now - 100
        else:
            lastfile = files[-1]
            lasthr = int(lastfile[-8:-3]) # most recent file
            nc_time_start = lasthr + 1
            numfiles = hours_now - lasthr
        

        print(nc_time_start)
        print(numfiles)
        numfiles = 1


        exec(open('eddy1_identification.py').read()) # does eddy analysis and debug plots
        exec(open('eddy2_track_eddies.py').read()) # tracks the eddys through time
        exec(open('eddy3_analyze_eddy_tracks.py').read()) # output of last timestep to csv for web graphics

        # uploads results to server
        session = ftplib.FTP('ftp_ip_add','ftp_user','ftp_pass') # ftp server
        file = open(fdir + file1,'rb')                  # file to send
        session.storbinary('STOR ' + file1, file)       # send the file
        file = open(fdir + file2,'rb')                 
        session.storbinary('STOR ' + file2, file)     
        file.close()                                    # close file and FTP
        session.quit()
        
        
        
    time.sleep(60*sleeptime) # sleep until next run

