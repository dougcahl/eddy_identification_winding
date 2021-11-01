#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 5 12:03:02 2021

@author: doug
"""
# call with master1.py avg_25hr radar_km_resolution region_num
# where 
# region_num = 1    region = 'USEGC'    # 1,2,6 km                              fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'
# region_num = 2    region = 'GAK'      # 1,2,6 km                              fullregionname = 'HFRADAR_US_Gulf_of_Alaska'
# region_num = 3    region = 'AKNS'     # 6 km                                  fullregionname = 'HFRADAR_Alaska_-_North_Slope'
# region_num = 4    region = 'USWC'     # 500m 1,2,6 km (500m is not supported)    fullregionname = 'HFRADAR_US_West_Coast'
# region_num = 5    region = 'USHI'     # 1,2,6 km                              fullregionname = 'HFRADAR_US_Hawaii'
# region_num = 6    region = 'PRVI'     # 2,6 km             fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands'
# example: python3 master1.py 1 6 1


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
from scipy.io import savemat
import re


np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)   



sleeptime = 10 # sleep between analysis in minutes
debug_plot = 0 # change to 0 when deployed

region_num = int(sys.argv[3])


while True:
#    for i_km in range(1,3):
#        if i_km == 0:
#            radar_km_resolution = 6         #% 1, 2 or 6 km nc file
#        if i_km == 1:
#            radar_km_resolution = 2         #% 1, 2 or 6 km nc file
#        if i_km == 2:
#            radar_km_resolution = 1         #% 1, 2 or 6 km nc file


#    avg_25hr                = 1         # set to 0 to use hourly currents. set to 1 for 25 hour average currents
    avg_25hr                = int(sys.argv[1])
#    radar_km_resolution     = 6         # 1, 2 or 6 km nc file, does not support west coast 500m
    radar_km_resolution     = int(sys.argv[2])
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
    # 1 and 2km grid sizes are converted below
    # grid_deg_size           = 2/3         # utm grid size for streamlines (ex. 1x1 grid is about 100km x 100km )
    # grid_deg_size           = 2/6         # utm grid size for streamlines (ex. 1x1 grid is about 100km x 100km )
    # scaling down grid size for 1 and 2 km data sets
    grid_deg_size = grid_deg_size*radar_km_resolution/6


    ### region selection
    if region_num == 1:
        region = 'USEGC'    # 1,2,6 km
        fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'
    if region_num == 2:
        region = 'GAK'      # 1,2,6 km
        fullregionname = 'HFRADAR_US_Gulf_of_Alaska'
    if region_num == 3:
        region = 'AKNS'     # 6 km
        fullregionname = 'HFRADAR_Alaska_-_North_Slope'
    if region_num == 4:
        region = 'USWC'     # 500m 1,2,6 km (500m is not supported)
        fullregionname = 'HFRADAR_US_West_Coast'
    if region_num == 5:
        region = 'USHI'     # 1,2,6 km
        fullregionname = 'HFRADAR_US_Hawaii'
    if region_num == 6:
        region = 'PRVI'     # 2,6 km
        fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands'



    # name setup
    if avg_25hr == 1:
        hrname = '25hr'
    else:
        hrname = 'hourly'


# done in bash script now
#    # create log directory if missing
#    dlog = 'logs/' # fiename in hours since 2019-09-29 20:00:00.000 UTC
#    logfn = region + str(radar_km_resolution) + hrname + '.txt'
#    if not os.path.isdir(dlog):
#        os.mkdir(dlog) # create directory


    # create Data directory if missing
    din = 'Data/' # fiename in hours since 2019-09-29 20:00:00.000 UTC
    if not os.path.isdir(din):
        os.mkdir(din) # create directory

    # directory where hf currents are saved
    din = 'Data/analysis_data/' # fiename in hours since 2019-09-29 20:00:00.000 UTC
    if not os.path.isdir(din):
        os.mkdir(din) # create directory
    din = din + region + '/'
    if not os.path.isdir(din):
        os.mkdir(din) # create directory
    din = din + hrname + '/'
    if not os.path.isdir(din):
        os.mkdir(din) # create directory
    din = din + 'km' + str(radar_km_resolution) + '/'
    if not os.path.isdir(din):
        os.mkdir(din) # create directory

    # results directory for ftp upload
    fdir = 'Data/eddy_tracks/'
    if not os.path.isdir(fdir):
        os.mkdir(fdir) # create directory
    fdir = fdir + region + '/'
    if not os.path.isdir(fdir):
        os.mkdir(fdir) # create directory
    fdir = fdir + hrname + '/'
    if not os.path.isdir(fdir):
        os.mkdir(fdir) # create directory
    fdir = fdir + 'km' + str(radar_km_resolution) + '/'
    if not os.path.isdir(fdir):
        os.mkdir(fdir) # create directory

#    file1 = region + '_eddy_tracks_' + str(radar_km_resolution) + 'km.csv'
#    file2 = region + '_eddy_tracks_' + str(radar_km_resolution) + 'km.npz'
#    file3 = region + '_eddy_tracks_' + str(radar_km_resolution) + 'km.mat'


    # get hourly offset and current hour of data
    if avg_25hr == 1:
        hr_name = '25hr'
        hr_name_full = '25_Hour_Average'
    else:
        hr_name = 'hourly'
        hr_name_full = 'Hourly'

    url1 = 'http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/' + region + \
    '/' + str(radar_km_resolution) + \
    'km/' + hr_name + \
    '/RTV/' + fullregionname + \
    '_' + str(radar_km_resolution) + \
     'km_Resolution_' + hr_name_full + \
    '_RTV_best.ncd.html'

    # get data set lat lon size    
    try:
        fp = urllib.request.urlopen(url1)
        mybytes = fp.read()
        mystr = mybytes.decode("utf8")
        fp.close()
    except OSError:
        print('website not giving description')
        continue

    # time offset
    # print(mystr)
    i_str = mystr.find('units: hours since ')
    mystr1_time = mystr[i_str+19:i_str+46];
#    print(mystr1_time) # offset from date
    # hours
    i_str = mystr.find('Float64 time_run[time =')
    mystr1 = mystr[i_str+24:i_str+32];
    i_str = mystr1.find('];')
    tnum = int(mystr1[0:i_str]) # latest time + 1, to set numfiles for this loop
#    nc_time = tnum - 1; 
    yyyy = int(mystr1_time[0:4])
    mm = int(mystr1_time[5:7])
    dd = int(mystr1_time[8:10])
    hh = int(mystr1_time[11:13])

    now = datetime.datetime.utcnow()
    start = datetime.datetime(yyyy,mm,dd)
    diff = now - start
    days, seconds = diff.days, diff.seconds
    hours_now = days * 24 + seconds // 3600
    hours_now = hours_now - hh; # change to hours_now = hours_now - 20-100;



#    now = datetime.datetime.utcnow()
#    start = datetime.datetime(2019,9,29)
#    diff = now - start
#    days, seconds = diff.days, diff.seconds
#    hours_now = days * 24 + seconds // 3600
#    hours_now = hours_now - 20; # change to hours_now = hours_now - 20-100;



    print('current time = ' + str(hours_now))
    if avg_25hr == 1:
        print('25 hr average 12.5 hours behind')


    files = sorted(os.listdir(din))
    if len(files) == 0:
        numfiles = 1
        nc_time_start = tnum - 20
        print('first run, starting 20 hours back from last data')
    else:
        lastfile = files[-1]
#        lasthr = int(lastfile[-8:-3]) # most recent file
        lasthr_str = re.split(r'[_.\s]\s*',lastfile)
        lasthr = int(lasthr_str[2])
        nc_time_start = lasthr + 1
        numfiles = hours_now - lasthr
        print('last file +1 = ' + str(nc_time_start))
        print('catchup num  = ' + str(numfiles))

#    numfiles = 1


    exec(open('eddy1_identification.py').read()) # does eddy analysis and debug plots
    exec(open('eddy2_track_eddies.py').read()) # tracks the eddys through time
    exec(open('eddy3_analyze_eddy_tracks.py').read()) # output of last timestep to csv for web graphics
#    break

    # uploads results to server
    # print('sending results')
    # try:
    #     session = ftplib.FTP('ftp_ip','ftpuser','ftp_pass') # ftp server
    #     file = open(fdir + file1,'rb')                  # file to send
    #     session.storbinary('STOR ' + file1, file)       # send the file
    #     file = open(fdir + file2,'rb')
    #     session.storbinary('STOR ' + file2, file)
    #     file = open(fdir + file3,'rb')
    #     session.storbinary('STOR ' + file3, file)
    #     file.close()                                    # close file and FTPsession = ftplib.FTP('104.128.235.148','ftp_home','homehome')
    #     session.quit()
    # except:
    #     print('ftp server down, no results going to webserver')


    print('sleep for ' + str(sleeptime) + ' minutes ...')
    time.sleep(60*sleeptime) # sleep until next run

