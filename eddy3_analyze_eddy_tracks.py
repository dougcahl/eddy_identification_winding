#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 12:03:02 2020

@author: doug
"""

#radar_km_resolution = 6;    # 1, 2 or 6 km nc file 

eddy_track_dist_param = 50; # 4.0; % distance in km for eddy to be continuation
eddy_track_time_param = 5; # how many hours can there be between identified eddies for it to be considered the same eddy


# directory root
din = 'Data/';

# create 'eddy_tracks' directory
din_data = din + 'eddy_tracks/'
din_data = din_data + 'km' + str(radar_km_resolution) + '/'



if not os.path.isdir(din_data):
    os.mkdir(din_data) # create directory


# ncout_name = din_data + 'eddy_tracks_' + str(radar_km_resolution) + 'km.nc'
npz_in_name = din_data + 'eddy_tracks_' + str(radar_km_resolution) + 'km.npz'
csv_out_name = din_data + 'eddy_tracks_' + str(radar_km_resolution) + 'km.csv'
csv_out_name_junk = din_data + 'eddy_tracks_' + str(radar_km_resolution) + 'km_junk.csv'




# a = np.load('/data/user/Files/Data/eddy_tracks/2020/eddy_tracks_6km.npz', allow_pickle=True)
a = np.load(npz_in_name, allow_pickle=True)



total_tracks           = a['total_tracks']
eddy_center_lon        = a['eddy_center_lon']
eddy_center_lat        = a['eddy_center_lat']
eddy_angular_vel       = a['eddy_angular_vel']
eddy_length_x          = a['eddy_length_x']
eddy_length_y          = a['eddy_length_y']
eddy_ellipse_theta     = a['eddy_ellipse_theta']
eddy_dir               = a['eddy_dir']
eddy_streamlines       = a['eddy_streamlines']
Time                   = a['Time']
eddy_track_dist_param  = a['eddy_track_dist_param']
eddy_track_time_param  = a['eddy_track_time_param']




time_all = Time.flatten
for tracki in range(0,int(total_tracks)):
    timetest = Time[tracki]
    if timetest.size > 1:
        timetest = timetest[-1]
    if tracki == 0:
        time_latest = timetest
    time_latest = np.maximum(timetest,time_latest)   





with open(csv_out_name_junk, 'w', newline='') as csvfileout:
    spamwriter = csv.writer(csvfileout, delimiter=',')
    spamwriter.writerow(['eddy_center_lon','eddy_center_lat','eddy_angular_vel',
             'eddy_length_x','eddy_length_y','eddy_ellipse_theta','eddy_dir','eddy_streamlines'])



with open(csv_out_name, 'w', newline='') as csvfileout,         open(csv_out_name_junk, 'r', newline='') as csvfilein:
    spamwriter = csv.writer(csvfileout, delimiter=',')
    row1 = csv.reader(csvfilein, delimiter=',')
    for row in row1:
        for i in range(0,100):
            rowlon ='lon' + str(i+1)
            row.append(rowlon)
            rowlat ='lat' + str(i+1)
            row.append(rowlat)
        spamwriter.writerow(row)




# now given a Time
# time1 = 737516.125
time1 = time_latest
for tracki in range(0,int(total_tracks)):
    # print(tracki)
    timetest = Time[tracki] # array of time values for this track
    if time1 in timetest: # if this time exists for this track
        timei = np.where(timetest == time1)
        timei = int(timei[0])
        # print('track ' + str(tracki) + '     timeindex ' + str(timei))

        if eddy_center_lon[tracki].size > 1: # eddy has a history or future from this timestep
            # positions 
            lons = eddy_center_lon[tracki][0:timei+1]
            lats = eddy_center_lat[tracki][0:timei+1]
            print(eddy_center_lon[tracki][0:timei+1]) #  -- past until this timestep
            print(eddy_center_lat[tracki][0:timei+1]) #  -- past until this timestep
            
            # eddy size for this timestep
            avel = eddy_angular_vel[tracki][timei]
            lx = eddy_length_x[tracki][timei]
            ly = eddy_length_y[tracki][timei]
            th = eddy_ellipse_theta[tracki][timei]
            direc = eddy_dir[tracki][timei]
            streamlines = eddy_streamlines[tracki][timei]
            
            print(eddy_angular_vel[tracki][timei])
            print(eddy_length_x[tracki][timei])
            print(eddy_length_y[tracki][timei])
            print(eddy_ellipse_theta[tracki][timei])
            print(eddy_dir[tracki][timei])
            print(eddy_streamlines[tracki][timei])
            
            
            with open(csv_out_name_junk, 'w', newline='') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=',')
                spamwriter.writerow([lons[0],lats[0],
                                    avel,lx,ly,th,direc,streamlines])
            
            
            with open(csv_out_name, 'a', newline='') as csvfileout,                     open(csv_out_name_junk, 'r', newline='') as csvfilein:
                spamwriter = csv.writer(csvfileout, delimiter=',')
                row1 = csv.reader(csvfilein, delimiter=',')
                for row in row1:
                    for i in range(0,lons.size-1):
                        rowlon = lons[i+1]
                        row.append(rowlon)
                        rowlat = lats[i+1]
                        row.append(rowlat)
                    spamwriter.writerow(row)
        
                    
            
            

        else: # eddy only exists at this timestep
            # positions 
            lons = eddy_center_lon[tracki]
            lats = eddy_center_lat[tracki]
            
            # eddy size for this timestep
            avel = eddy_angular_vel[tracki]
            lx = eddy_length_x[tracki]
            ly = eddy_length_y[tracki]
            th = eddy_ellipse_theta[tracki]
            direc = eddy_dir[tracki]
            streamlines = eddy_streamlines[tracki]
            
            with open(csv_out_name, 'a', newline='') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=',')
                spamwriter.writerow([lons,lats,
                                    avel,lx,ly,th,direc,streamlines])
            
       
            
            
            # positions 
            print(eddy_center_lon[tracki])
            print(eddy_center_lat[tracki])
            
            # eddy size for this timestep
            print(eddy_angular_vel[tracki])
            print(eddy_length_x[tracki])
            print(eddy_length_y[tracki])
            print(eddy_ellipse_theta[tracki])
            print(eddy_dir[tracki])
            print(eddy_streamlines[tracki])

    else:
        print('this time not in this track')

