#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 12:03:02 2020

@author: doug
"""
#radar_km_resolution = 6;    # 1, 2 or 6 km nc file 

#eddy_track_dist_param = 25; # 4.0; % distance in km for eddy to be continuation
#eddy_track_time_param = 5; # how many hours can there be between identified eddies for it to be considered the same eddy


# directory root
din = 'Data/';

# import data directory
din_data = din + 'analysis_data/'
din_data = din_data + region + '/'
din_data = din_data + hrname + '/'
din_data = din_data + 'km' + str(radar_km_resolution) + '/'

# create 'eddy_tracks' directory
dout_data = din + 'eddy_tracks/'
if not os.path.isdir(dout_data):
    os.mkdir(dout_data) # create directory

# save data directory for year
dout_data = dout_data + region + '/'
if not os.path.isdir(dout_data):
    os.mkdir(dout_data) # create directory

dout_data = dout_data + hrname + '/'
if not os.path.isdir(dout_data):
    os.mkdir(dout_data) # create directory

dout_data = dout_data + 'km' + str(radar_km_resolution) + '/'
if not os.path.isdir(dout_data):
    os.mkdir(dout_data) # create directory

# ncout_name = dout_data + 'eddy_tracks_' + str(radar_km_resolution) + 'km.nc'
npz_out_name = dout_data + region + '_eddy_tracks_' + str(radar_km_resolution) + 'km.npz'
mat_out_name = dout_data + region + '_eddy_tracks_' + str(radar_km_resolution) + 'km.mat'

# identified eddies for it to be considered the same eddy
files = sorted(os.listdir(din_data)) # get filenames
file_count = 0

# initizalize vars
eddy_center_lat     = []
eddy_center_lon     = []
eddy_dir            = []
eddy_angular_vel    = []
eddy_streamlines    = []
eddy_length_x       = []
eddy_length_y       = []
eddy_ellipse_theta  = []
ttime                = []
eddy_num            = [] # number of eddies for each timestep




total_eddies = 0 # total eddies
for x in files:
    # x = '201701041800_hfr_usegc_6km_rtv_uwls_NDBC.nc'
    if x.find(str(radar_km_resolution) + 'km') < 0: # resolution selection
        continue
    
#    print(x) # print filename
    
    # read nc data
    ncin_name = din_data + x
    ncin = Dataset(ncin_name, mode='r') # open nc file
    
    
    # get data
    eddy_center_lat.append(np.array(ncin.variables['eddy_center_lat'][:]))
    eddy_center_lon.append(np.array(ncin.variables['eddy_center_lon'][:]))
    eddy_dir.append(np.array(ncin.variables['eddy_dir'][:]))
    eddy_angular_vel.append(np.array(ncin.variables['eddy_angular_vel'][:]))
    eddy_streamlines.append(np.array(ncin.variables['eddy_streamlines'][:]))
    eddy_length_x.append(np.array(ncin.variables['eddy_length_x'][:]))
    eddy_length_y.append(np.array(ncin.variables['eddy_length_y'][:]))
    eddy_ellipse_theta.append(np.array(ncin.variables['eddy_ellipse_theta'][:]))
    eddies = len(np.array(ncin.variables['eddy_dir'][:]))
    ncin.close() # close nc file
    
    
    # timestamp of file
    # yyyy,    mm,    dd,     HH,    MM
    # x[0:4],x[4:6],x[6:8],x[8:10],x[10:12]
    
    # hours since 2019-09-29 20:00:00.000 UTC
    # UTC time, add 366 to get matlab time
    ttime.append(datetime.datetime(2019, 9, 29, 20, 0) + datetime.timedelta(hours = int(x[-8:-3])))
    #time.append(date.toordinal(date(int(x[0:4]),int(x[4:6]),int(x[6:8]))) + float(x[8:10])/24 + float(x[10:12])/24/60) # add 366 to get matlab time
    
    eddy_num.append(eddies) # eddies at this timestep
    total_eddies = total_eddies + eddies # total eddies found


#print(eddy_center_lon)
#print(eddy_center_lat)
# if not on first timestep
# check each eddy on this timestep for distance from previous eddies
# if distance < eddy_track_dist_param
# andif same winding dir as before
# then this eddy is a continuation of the last eddy
# check if this eddy track exists
# if it does, add new data point
# else look back eddy_track_time_param hours/timesteps for eddy
# else if no previous eddy, make new eddy track


print('Total eddies = ' + str(total_eddies))
print('Max eddies in a single timestep = ' + str(np.max(eddy_num)))

# no tracking for first timestep
it = 0; # timestep 1
i_cluster = eddy_num[it]

# initialize vars since we're not in Matlab anymore
lat_center      = []
lon_center      = []
ellipse_theta   = []
eig1            = []
eig2            = []
omega           = [] # rads per sec
direction       = [] # +1 is cyclonic, -1 is anticyclonic
num_streams     = [] # intensity, number of streamliens
Time            = []
timegap         = []

# loop through eddies found in this timestep
for i in range(0,i_cluster):
    # print(i)
    lat_center.append(np.array(eddy_center_lat[it][i]))
    lon_center.append(np.array(eddy_center_lon[it][i]))
    ellipse_theta.append(np.array(eddy_ellipse_theta[it][i]))
    eig1.append(np.array(eddy_length_x[it][i]))
    eig2.append(np.array(eddy_length_y[it][i]))
    omega.append(np.array(eddy_angular_vel[it][i])) # rads per sec
    direction.append(np.array(eddy_dir[it][i])) # +1 is cyclonic, -1 is anticyclonic
    num_streams.append(np.array(eddy_streamlines[it][i])) # intensity, number of streamliens
    Time.append(np.array(ttime[it]));
    timegap.append(np.array(1));


# 
total_tracks = i_cluster;
#print('Starting track finding routine')
#print('Total number of timestep to analyze = ' + str(len(ttime)))
#print('Timestep = ' + str(it))

# track the same eddies
for it in range(1,len(ttime)):
#    print('Timestep = ' + str(it))
    i_cluster = len(eddy_center_lat[it]); # number of eddies at this timestep
    for i in range(0,i_cluster):
        curr_lat_center = eddy_center_lat[it][i];
        curr_lon_center = eddy_center_lon[it][i];
#        print(curr_lat_center)
#        print(curr_lon_center)
        #         currentpos = [x(it,i) y(it,i)];
        #         cur_size = [eigvals1(it,i) eigvals2(it,i)];
        # check over all tracks with timesteps from before
        succ = 0
        # i1_start = find(timegap <= eddy_track_time_param,1);
        # potential_tracks = np.argwhere(np.array(timegap) <= eddy_track_time_param)
        
        # timegaps less than eddy_track_time_param and >0 to avoid eddies in this timestep from merging
        potential_tracks = np.argwhere((np.array(timegap) <= eddy_track_time_param) & (np.array(timegap) > 0))
        for ii1 in potential_tracks:
            i1 = int(ii1)
            # calculate distance between eddy centers in km
            if np.size(lat_center[i1]) == 1:
                dx,dy,zone_a,zona_a_str = utm.from_latlon(np.array([curr_lat_center, lat_center[i1]]),
                                                          np.array([curr_lon_center, lon_center[i1]]))
                dist = np.sqrt((dx[1]-dx[0])**2 + (dy[1]-dy[0])**2)/1000 # in km
                
            else:
                dx,dy,zone_a,zona_a_str = utm.from_latlon(np.array([curr_lat_center, lat_center[i1][-1]]),
                                                          np.array([curr_lon_center, lon_center[i1][-1]]))
                dist = np.sqrt((dx[1]-dx[0])**2 + (dy[1]-dy[0])**2)/1000 # in km
            
          
            if dist < eddy_track_dist_param:
                if timegap[i1] <= eddy_track_time_param:
                    lat_center[i1]    = np.append(lat_center[i1] , eddy_center_lat[it][i])
                    lon_center[i1]    = np.append(lon_center[i1] , eddy_center_lon[it][i])
                    ellipse_theta[i1]       = np.append(ellipse_theta[i1] , eddy_ellipse_theta[it][i])
                    eig1[i1]        = np.append(eig1[i1] , eddy_length_x[it][i])
                    eig2[i1]        = np.append(eig2[i1] , eddy_length_y[it][i])
                    omega[i1]       = np.append(omega[i1] , eddy_angular_vel[it][i])
                    direction[i1]   = np.append(direction[i1] , eddy_dir[it][i])  # +1 is cyclonic, -1 is anticyclonic
                    num_streams[i1] = np.append(num_streams[i1] , eddy_streamlines[it][i])  # intensity
                    Time[i1]        = np.append(Time[i1] , ttime[it])
#                    print('Time gaps away = ' + str(timegap[i1]))
                    timegap[i1]     = 0; # zero because it's gonna add one right after this loop
                    succ = 1;
                    break; # do not combine more than 2 tracks
                
 
        # if not close to others, create new track
        if succ == 0:
            total_tracks = total_tracks + 1;
            lat_center.append(np.array(eddy_center_lat[it][i]));
            lon_center.append(np.array(eddy_center_lon[it][i]));
            ellipse_theta.append(np.array(eddy_ellipse_theta[it][i]));
            eig1.append(np.array(eddy_length_x[it][i]));
            eig2.append(np.array(eddy_length_y[it][i]));
            omega.append(np.array(eddy_angular_vel[it][i]));
            direction.append(np.array(eddy_dir[it][i])); # +1 is cyclonic, -1 is anticyclonic (in Northern Hemisphere)
            num_streams.append(np.array(eddy_streamlines[it][i])); # intensity
            Time.append(np.array(ttime[it]));
            timegap.append(np.array(0)); # zero because it's gonna add one right after this loop

    
    # add 1 timestep to the counter
    for i1 in range(0,total_tracks):
        timegap[i1] = timegap[i1] + 1;
    


         

# save nc
# convert to numpy arrays
total_tracks           = np.array(total_tracks,dtype=object) # eddy center lon/lat
eddy_center_lat        = np.array(lat_center,dtype=object)
eddy_center_lon        = np.array(lon_center,dtype=object) # eddy direction, math notation
eddy_angular_vel       = np.array(omega,dtype=object) # angular velocty (deg/sec) math notation of eddy
eddy_length_x          = np.array(eig1,dtype=object) # number of streamlines for each eddy
eddy_length_y          = np.array(eig2,dtype=object) # ellipse fit x length (meters)
eddy_ellipse_theta     = np.array(ellipse_theta,dtype=object) # ellipse fit y length (meters)
eddy_dir               = np.array(direction,dtype=object) # rotation of eddy ellipse fit, math degrees
eddy_streamlines       = np.array(num_streams,dtype=object) 
Time                   = np.array(Time,dtype=object) 
# timegap = np.array(timegap) 
eddy_track_dist_param  = np.array(eddy_track_dist_param,dtype=object) 
eddy_track_time_param  = np.array(eddy_track_time_param,dtype=object) 

Time_mat = []
for i in range(0,Time.size):
    Time_1 = Time[i]
    Time_mat2 = []
    for j in range(0,Time_1.size):
        if Time_1.size == 1:
            Time_2 = Time_1.tolist()
        else:
            Time_2 = Time_1[j]
        Time_mat2.append(Time_2.strftime("%d-%b-%Y %H:%M:%S"))
        #        Time_mat[i][j] = Time_2.strftime("%d-%b-%Y %H:%M:%S")
    Time_mat.append(Time_mat2)


    

np.savez(npz_out_name, 
         total_tracks          = total_tracks,
         eddy_center_lat       = eddy_center_lat,
         eddy_center_lon       = eddy_center_lon,
         eddy_angular_vel      = eddy_angular_vel,
         eddy_length_x         = eddy_length_x,
         eddy_length_y         = eddy_length_y,
         eddy_ellipse_theta    = eddy_ellipse_theta,
         eddy_dir              = eddy_dir,
         eddy_streamlines      = eddy_streamlines,
         Time                  = Time,
         eddy_track_dist_param = eddy_track_dist_param,
         eddy_track_time_param = eddy_track_time_param)

# a = np.load('/home/doug/Documents/MATLAB/eddy/data/eddy_tracks/2017/eddy_tracks_6km.npz', allow_pickle=True)




mdic = {"total_tracks": total_tracks, 
         "eddy_center_lat": eddy_center_lat,
         "eddy_center_lon": eddy_center_lon,
         "eddy_angular_vel": eddy_angular_vel,
         "eddy_length_x": eddy_length_x,
         "eddy_length_y": eddy_length_y,
         "eddy_ellipse_theta": eddy_ellipse_theta,
         "eddy_dir": eddy_dir,
         "eddy_streamlines": eddy_streamlines,
         "Time": Time_mat,
         "eddy_track_dist_param": eddy_track_dist_param,
         "eddy_track_time_param": eddy_track_time_param }

savemat(mat_out_name, mdic) # save .mat file




