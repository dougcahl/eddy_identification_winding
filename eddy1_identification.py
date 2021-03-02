#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 12:03:02 2020

@author: doug
"""




din = 'Data/' # '2017/' # directory of netcdf files to process/save
# dout = '/home/doug/Documents/MATLAB/eddy/data/eddy_analysis_py/'
# debug_plot  = 1 # % debug plots
#debug_info  = 1 # print when eddy is found for each streamline
# numfiles    = 1 # % how many files to run. Set to 0 to run all files
#
#radar_km_resolution = 6         #% 1, 2 or 6 km nc file
#res_x_mult          = 1         #% for higher resolution, if >1 results in this squared processing time (i.e. 3 => 9 times longer), only implemented in Matlab
#winding_thres       = 270       #% how much winding (in degrees) to be a eddy
#d_thres             = 25        #% streamline starting and ending max distance in km
#baddir_thres        = 15        #% breaks at 15 degrees of bad dirs
#param_center_dist   = 25        #% how many km are clusters together
#min_pts             = 5         # min number of grid points streamline has
#
# % utm grid size for streamlines (ex. 1x1 grid is about 100km x 100km )
#grid_deg_size = 2;
#
#region = 'USEGC'
# nc_time = 11764 # hours since 2019-09-29 20:00:00.000 UTC
# nc_time_start = 11800 # 2021-1-31 21:00:00 UTC
# nc_time_start = 12000 # 2021-1-31 21:00:00 UTC
# nc_time_start = 12303 # 2021-1-31 21:00:00 UTC
#nc_time_start = 11801 # 2021 ??? UTC
# time.append(date.toordinal(date(int(x[0:4]),int(x[4:6]),int(x[6:8]))) + float(x[8:10])/24 + float(x[10:12])/24/60) # add 366 to get matlab time
#time.append(date.toordinal(date(int(x[0:4]),int(x[4:6]),int(x[6:8]))) + float(x[8:10])/24 + float(x[10:12])/24/60) # add 366 to get matlab time
#start_time_ordinal = date.toordinal(date(2019,9,29)) + 20/24 + nc_time_start/24 # UTC time, add 366 to get matlab time
start_time_ordinal = datetime.datetime(2019, 9, 29, 20, 0) + datetime.timedelta(hours = nc_time_start)




# create output directory structure
# if not os.path.isdir(dout):
    # os.mkdir(dout) # create directory

# output directory for saving analysis data
dout_data = din + 'analysis_data/'
if not os.path.isdir(dout_data):
    os.mkdir(dout_data) # create directory

dout_data = dout_data + 'km' + str(radar_km_resolution)  + '/'
if not os.path.isdir(dout_data):
    os.mkdir(dout_data) # create directory
    
# # output directory for saving analysis data
# dout_data = din + 'analysis_data/' + str(yr) + '/'
# if not os.path.isdir(dout_data):
#     os.mkdir(dout_data) # create directory

# output directory for saving analysis figures
dout_fig = din + 'analysis_figs/'
if not os.path.isdir(dout_fig):
    os.mkdir(dout_fig) # create directory

dout_fig = dout_fig + 'km' + str(radar_km_resolution) + '/'
if not os.path.isdir(dout_fig):
    os.mkdir(dout_fig) # create directory
    
# # output directory for saving analysis figures
# dout_fig = din + 'analysis_figs/' + str(yr) + '/'
# if not os.path.isdir(dout_fig):
#     os.mkdir(dout_fig) # create directory

# # input directory
# din = din + str(yr) + '/'
# files = sorted(os.listdir(din)) # get filenames





# start processing
file_count = 0
while file_count < numfiles:
    nc_time = nc_time_start + file_count # hours since 2019-09-29 20:00:00.000 UTC
    # nc_time_ordinal = date.toordinal(date(2019,9,29)) + 20/24 + nc_time_start/24) # UTC time, add 366 to get matlab time
    nc_time_ordinal = datetime.datetime(2019, 9, 29, 20, 0) + datetime.timedelta(hours = nc_time) # UTC time
# for x in files:
    # x = '201701010100_hfr_usegc_6km_rtv_uwls_NDBC.nc'
    # x = '201701012200_hfr_usegc_6km_rtv_uwls_NDBC.nc'
    
    
    
    # if x.find(str(radar_km_resolution) + 'km') < 0: # resolution selection
        # continue
    
    timer1 = time.time()
    file_count = file_count + 1
   
    # print(x) # print filename
    # ncin_name = din + x
    # ncout_name = dout_data + x[0:-3] + '.nc'
    # fig_name = dout_fig + x[0:-3]
    # ncin = Dataset(ncin_name, mode='r')
    
    print('running analysis on file ' + str(nc_time))
    fig_name = dout_fig + region + '_' + str(radar_km_resolution) + 'km_' + str(nc_time)
    
#     # testing 25 hr opendap dataset
#     url = 'http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USEGC/6km/25hr/RTV/HFRADAR_US_East_and_Gulf_Coast_6km_Resolution_25_Hour_Average_RTV_best.ncd?lat[0:1:459],lon[0:1:700],time[11746],u_mean[11746][0:1:459][0:1:700],v_mean[11746][0:1:459][0:1:700]'
    # region = 'USEGC'
    # nc_time = 11764 # hours since 2019-09-29 20:00:00.000 UTC
    ncout_name = dout_data + region + '_' + str(radar_km_resolution) + 'km_' + str(nc_time) + '.nc' 
    url = 'http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/' + region + '/' + str(radar_km_resolution) +         'km/25hr/RTV/HFRADAR_US_East_and_Gulf_Coast_' + str(radar_km_resolution) + 'km_Resolution_25_Hour_Average_RTV_best.ncd?' +         'lat[0:1:459],lon[0:1:700],time[' + str(nc_time) + '],u_mean[' + str(nc_time) +         '][0:1:459][0:1:700],v_mean[' + str(nc_time) + '][0:1:459][0:1:700]'
    
    try:
        print('downloading data ...')
        ncin = Dataset(url)
    except OSError:
        print('file does not exist or website down')
        print(url)
        break
    
    print('running eddy analysis ...')
    lon = ncin.variables['lon'][:]
    lat = ncin.variables['lat'][:]
    t = ncin.variables['time'][:] # hours since 2019-09-29 20:00:00.000 UTC
    u = ncin.variables['u_mean'][:]
    v = ncin.variables['v_mean'][:]
    u.mask = np.isnan(u)

    # standard read files
#     lon = ncin.variables['lon'][:]
#     lat = ncin.variables['lat'][:]
#     t = ncin.variables['time'][:] # matlab time?
#     u = ncin.variables['u'][:]
#     v = ncin.variables['v'][:]
    
    ncin.close()
    
    

    
    # print(t)
    lons,lats = np.meshgrid(lon, lat)
    u1 = np.squeeze(u)                     
    v1 = np.squeeze(v)                     

    
#     # quiver plot of data, uses a lot of ram 
#     if debug_plot == 1:
#         # plt.figure()
#         # plt.plot(lon) # 20:40
#         # plt.plot(lat) # -100:-60
        
#         fig = plt.figure(figsize=(20,12))
#         # plt.quiver(lons,lats,u1,v1)
#         plt.quiver(lons[~u1.mask],lats[~u1.mask],u1[~u1.mask],v1[~u1.mask])
#         plt.axis([np.min(lon),np.max(lon),np.min(lat),np.max(lat)])
#         plt.pause(0.05)
#         plt.savefig(fig_name + '_wholemap.png')
#         # plt.axis('equal')
#         # fig.canvas.draw()
#         # fig.canvas.flush_events()
    
    
    # initialize these to fill later
    i_streams = -1
    streamsx = [] # np.array()
    streamsy = [] # np.array()
    streamsx_center = []
    streamsy_center = []
    winding_dir = []
    angular_vel = []
    streams_lat = []
    streams_lon = []
    streams_center_lat = []
    streams_center_lon = []
    streams_ZONE = []
    streams_ZONE_str = []
                    
    # start calculating streamlines for this file
    count = 1
    # break1 = False
    for i in range(0,u1.shape[0]): # 1 - 460
        for j in range(0,u1.shape[1]): # 1 - 701
            if not u1.mask[i,j]: # if the data is not masked, run streamline
                # print(count)
                count = count + 1
                
                # start position of streamline
                start_lon = lons[i,j]
                start_lat = lats[i,j]
                
                # get an area grid_deg_size around it
                startx_high = start_lon + grid_deg_size/2
                startx_low  = start_lon - grid_deg_size/2
                starty_high = start_lat + grid_deg_size/2
                starty_low  = start_lat - grid_deg_size/2
                
                
                i_grid1 = lons >= startx_low
                i_grid2 = lons <= startx_high
                i_grid3 = lats >= starty_low
                i_grid4 = lats <= starty_high
                i_grid12 = np.logical_and(i_grid1, i_grid2)
                i_grid34 = np.logical_and(i_grid3, i_grid4)
                i_grid = np.logical_and(i_grid12, i_grid34)
                  
                lons_grid = lons[i_grid]
                lats_grid = lats[i_grid]
                u_grid = u1[i_grid]
                v_grid = v1[i_grid]
                
                # check start is same zone
                # startx, starty, ZONE_start, ZONE_start_str = utm.from_latlon(start_lat, start_lon)
                x_grid,y_grid, ZONE_start, ZONE_start_str = utm.from_latlon(lats_grid,lons_grid)
                startx = np.mean(x_grid);
                starty = np.mean(y_grid);
                
                d_grid = radar_km_resolution*1000/res_x_mult
                xgrid = np.arange(np.min(x_grid),np.max(x_grid)+d_grid/2,d_grid)
                ygrid = np.arange(np.min(y_grid),np.max(y_grid)+d_grid/2,d_grid)
                xgrid,ygrid = np.meshgrid(xgrid,ygrid)
                
                
               
                
                # interpolate onto equally spaced grid, is this necessary for python streamline?
                # points = np.array([x_grid[~u_grid.mask],y_grid[~u_grid.mask]])
                # points = np.transpose(points)
                # ugrid1 = griddata(points, u_grid[~u_grid.mask], (xgrid.flatten(),ygrid.flatten()), method='linear')
                # vgrid1 = griddata(points, v_grid[~u_grid.mask], (xgrid.flatten(),ygrid.flatten()), method='linear')
                # ugrid11 = np.reshape(ugrid1,xgrid.shape)
                # vgrid11 = np.reshape(vgrid1,xgrid.shape)
                
                # interpolate onto equally spaced grid
                u_grid_0 = u_grid;
                v_grid_0 = v_grid;
                u_grid_0[u_grid.mask] = np.nan
                v_grid_0[u_grid.mask] = np.nan
                points = np.array([x_grid,y_grid])
                points = np.transpose(points)
                ugrid1 = griddata(points, u_grid_0, (xgrid.flatten(),ygrid.flatten()), method='linear')
                vgrid1 = griddata(points, v_grid_0, (xgrid.flatten(),ygrid.flatten()), method='linear')
                ugrid11 = np.reshape(ugrid1,xgrid.shape)
                vgrid11 = np.reshape(vgrid1,xgrid.shape)
                
                # calculate streamline
                points = streamplot_dc(xgrid,ygrid,ugrid11,vgrid11,
                               linewidth=2,
                               cmap=plt.cm.coolwarm,start_points=np.array([[startx],[starty]]).T,
                               integration_direction='forward')
                
                if len(points) < 10: # if no streamline
                    continue
                
                # else do winding angle for this streamline
                stline_x = points[:,0,0] # x points of streamline
                stline_y = points[:,0,1] # y points of streamline
                 
                # initial angle
                eddy        = 0 #; % not an eddy initially
                winding     = 0 #; % winding angle
                dir1        = 0 #; % direction of winding
                baddir      = 0 #; % initially no bad directions
                
                
                # remove duplicate values
                st_len = len(stline_x)
                ii = 0
                while ii < st_len-1:
                    # print(ii)
                    if stline_x[ii+1] == stline_x[ii] and stline_y[ii+1] == stline_y[ii]:
                        # del stline_x[ii+1]
                        stline_x = np.concatenate((stline_x[:ii+1],stline_x[ii+2:]))
                        stline_y = np.concatenate((stline_y[:ii+1],stline_y[ii+2:]))
                        st_len = st_len - 1
                    
                    ii = ii + 1
                
                # first winding angle
                dx          = stline_x[1] - stline_x[0] #;
                dy          = stline_y[1] - stline_y[0] #;
                ang0        = np.arctan2(dy,dx)*180/math.pi #;
                
                # follow streamline
                for ii in range(1,len(stline_x)-1): # % no angle at last data point
                    dx = stline_x[ii+1] - stline_x[ii];
                    dy = stline_y[ii+1] - stline_y[ii];
                    ang1 = np.arctan2(dy,dx)*180/math.pi #; % new angle
                    ang_diff = ang1 - ang0 #; % winding angle is extra turning
                    ang0 = ang1 #; % reset angle
                    if ang_diff > 180:
                        ang_diff = ang_diff - 360 #;
                    if ang_diff < -180:
                        ang_diff = ang_diff + 360 #;
                    
                    winding = winding + ang_diff # sum the winding angle right away
                    # print(ang_diff)
                    if dir1 == 0: # % set with first point
                        dir1 = np.sign(ang_diff) #; # % if returns zero will set next time
                    
                    if ii > 1: # % only after first winding angle
                        newdir = np.sign(ang_diff) #;
#                         winding = winding + ang_diff # sum the winding angle after first angle is set (ii>1)
#                         newdir = np.sign(ang_diff) #;
                        if newdir != dir1: # % if winding changes direction
                            baddir = baddir + ang_diff #;
                            if np.abs(baddir) > baddir_thres: # % break after this many degrees bad dirs
                                break
                        if newdir == dir1: # % reset if goes back to correct direction before baddir_thres
                            baddir = 0 #;
                        
                    
                    if np.abs(winding) > winding_thres: # % if winding angle is close to 360
                        # % check distance from starting point
                        dx = stline_x[ii+1] - stline_x[0] #;
                        dy = stline_y[ii+1] - stline_y[0] #;
                        d = np.sqrt(dx**2+dy**2) #;
                        if eddy == 1: # % if already eddy see if distance is closer
                            if d > d1: # % if distance increases break and save
                                break
                            
                        
                        if d < d_thres*1000: # % now it is an eddy!
                            eddy = 1 #;
                            d1 = d #; % save this distance
                        # % we should still keep going to see if d decreases
                    # % otherwise keep going
                # for streamline
                
                
                
                
                
                
                # streamsx = np.array()
                # streamsy = np.array()
                # if it is an eddy
                # eddy = 1
                if eddy == 1 and ii > min_pts - 2:
                    i_streams = i_streams + 1;
                    # print(ii)
#                    if debug_info == 1:
#                        if i_streams < 1:
#                            print('found ' + str(i_streams + 1) + ' winding streamline!')
#                        else:
#                            print('found ' + str(i_streams + 1) + ' winding streamlines!')
                 
                    # streamline points
                    streamsx.append(stline_x[0:ii+1]); # % up to last pt before break
                    streamsy.append(stline_y[0:ii+1]);
                    
                    # streamline center
                    streamsx_center.append(np.mean(stline_x[0:ii+1]));
                    streamsy_center.append(np.mean(stline_y[0:ii+1]));
#                    print(stline_x[0:ii+1])
#                    print(stline_y[0:ii+1])
#                    print(np.mean(stline_x[0:ii+1]))
#                    print(np.mean(stline_y[0:ii+1]))
                    
                    # winding direction
                    winding_dir.append(dir1); # % math notation
                    
                    # % angular velocity in deg/sec
                    u1a = griddata((xgrid.flatten(),ygrid.flatten()), ugrid1,
                                   (streamsx[i_streams],streamsy[i_streams]), method='linear')
                    v1a = griddata((xgrid.flatten(),ygrid.flatten()), vgrid1,
                                   (streamsx[i_streams],streamsy[i_streams]), method='linear')
                    vv = np.sqrt(u1a**2 + v1a**2); # % in m/s
                    vv = vv[0:-1]; # % last velocity doesn't add to streamline
                    dx = np.diff(streamsx[i_streams]);
                    dy = np.diff(streamsy[i_streams]);
                    dr = np.sqrt(dx**2 + dy**2); # % in m
                    dtime = dr/vv; # % time in seconds
                    dtime = np.sum(dtime);
                    angular_vel.append(winding/dtime); # % degs/second
                    
                    # % convert to lat lon for streams
                    [streams_lat1,streams_lon1] = utm.to_latlon(streamsx[i_streams],streamsy[i_streams],
                                                                ZONE_start,ZONE_start_str);
                    streams_lat.append(streams_lat1)
                    streams_lon.append(streams_lon1)
                    
                    [streams_center_lat1,streams_center_lon1] = utm.to_latlon(streamsx_center[i_streams],
                                                  streamsy_center[i_streams],ZONE_start,ZONE_start_str);
                    streams_center_lat.append(streams_center_lat1)
                    streams_center_lon.append(streams_center_lon1)
#                    print(streams_center_lon1)
#                    print(streams_center_lat1)
                    
                    streams_ZONE.append(ZONE_start);
                    streams_ZONE_str.append(ZONE_start_str);

                   
                
                
                
               
                if debug_plot == 1 and eddy == 1:
                    # plot grid and equally spaced grid
                    # plt.figure()
                    # plt.plot(xgrid,ygrid,'ks')
                    # plt.plot(xgrid,ygrid,'r.')
                    # plt.plot(startx,starty,'bs')
                    
                    # utm quiver
                    # plt.figure()
                    # plt.quiver(x_grid,y_grid,u_grid,v_grid)
                    # plt.axis('equal')
                    
                    # quiver interpolated onto equally spaced grid
                    # plt.figure()
                    # plt.quiver(xgrid,ygrid,ugrid1,vgrid1)
                    # plt.axis('equal')
                    
                    # streamline on equally spaced grid
                    # plt.figure()
                    # plt.streamplot(xgrid,ygrid,ugrid11,vgrid11,
                    #                linewidth=2,
                    #                cmap=plt.cm.coolwarm)
                    
                    # # quiver plot with streamline on raw quivers
                    fig = plt.figure(figsize=(10,6))
                    plt.quiver(x_grid,y_grid,u_grid,v_grid)
                    plt.plot(stline_x[0:ii+1],stline_y[0:ii+1],'-bs')
                    plt.axis('equal')
                    plt.savefig(fig_name + '_streamline_' + str(i_streams) + '.png')
                    plt.close(fig)
                    
                    fig = plt.figure(figsize=(10,6))
                    plt.quiver(x_grid,y_grid,u_grid,v_grid)
                    plt.plot(stline_x[0:ii+1],stline_y[0:ii+1],'-bs')
                    plt.axis('equal')
                    plt.axis([np.min(stline_x[0:ii+1])-1000,np.max(stline_x[0:ii+1])+1000,
                              np.min(stline_y[0:ii+1])-1000,np.max(stline_y[0:ii+1])+1000])
                    plt.savefig(fig_name + '_streamline_zoomed' + str(i_streams) + '.png')
                    plt.close(fig)
                    # fig.canvas.draw()
                    # fig.canvas.flush_events()
                    # plt.draw()
                    
                    
                    # quiver plot with streamline, equispaced meshgrid
                    # plt.figure(figsize=(12, 8))
                    # plt.quiver(xgrid,ygrid,ugrid1,vgrid1)
                    # plt.plot(stline_x,stline_y,'-bs')
                    # plt.axis('equal')
                    
                    # plot of streamline
                    # plt.figure()
                    # plt.plot(stline_x,stline_y,'-bs')
                    # plt.axis('equal')
                    # break
                    # break1 = True
                    # break
           
        
        # break if above loop breaks
        # if break1:
        #     break
                
                
   
 
    
   # %% streamline clustering
    
    streams_center_lat = np.array(streams_center_lat)
    streams_center_lon = np.array(streams_center_lon)
    # streams_lat = np.array(streams_lat)
    # streams_lon = np.array(streams_lon)
    winding_dir = np.array(winding_dir)
    angular_vel = np.array(angular_vel)
    # streams_ZONE
    # streams_ZONE_str
    
    # eddy_nums = 1 :i_streams; # % eddy numbers
    eddy_nums = np.array(range(0,i_streams+1))
    for i1 in range(0,i_streams):
        if eddy_nums[i1] != i1:
            continue
        
        # % check to see if centers are close enough
        for j1 in range(i1,i_streams+1):
            if i1 == j1:
                continue
            
            dx,dy,zone_a,zona_a_str = utm.from_latlon(streams_center_lat[[i1,j1]],streams_center_lon[[i1,j1]])
            d = np.sqrt((dx[1]-dx[0])**2 + (dy[1]-dy[0])**2)/1000 # in km
            
            # % if centers are close and winding direction is the same
            if d < param_center_dist and (winding_dir[i1] == winding_dir[j1]): # % then eddy are in the same cluster
                eddy_nums[j1] = i1; # % same eddy as i
                # if debug_info == 1:
                    # print('same eddy') # % debug console info
            
        
    # break
    # %% save clusters here
    eddy_center_x       = []
    eddy_center_y       = []
    eddy_center_lat     = []
    eddy_center_lon     = []
    eddy_dir            = []
    eddy_angular_vel    = []
    eddy_streamlines    = []
    eddy_ellipse_theta  = []
    eddy_length_x       = []
    eddy_length_y       = []
    
    
    
    i1 = np.unique(eddy_nums);
    eddies = 0;
    # for j1 = 1:length(i1)
    for j1 in range(0,len(i1)):
        k = np.argwhere(eddy_nums == i1[j1]);
        eddies = eddies + 1;
        ptsx = [];
        ptsy = [];
        eddy_angular_vel1 = [];
        
        zone_eddy = streams_ZONE[k[0][0]]; # % UTM zone for this cluster
        zone_eddy_str = streams_ZONE_str[k[0][0]]; # % UTM zone for this cluster
        # for ii = 1:length(k)
        for ii in range(0,len(k)):
            eddy_angular_vel1.append(angular_vel[k[ii][0]])
            streamsx,streamsy,zone_a,zona_a_str = utm.from_latlon(streams_lat[k[ii][0]],
                      streams_lon[k[ii][0]],zone_eddy,zone_eddy_str)
            if len(ptsx) == 0:
                ptsx = streamsx
                ptsy = streamsy
            else:
                ptsx = np.concatenate((ptsx,streamsx))
                ptsy = np.concatenate((ptsy,streamsy))
            
#        print(zone_eddy)
#        print(zone_eddy_str)
#        print(zone_a)
#        print(zona_a_str)
        eddy_center_x.append(np.mean(ptsx))
        eddy_center_y.append(np.mean(ptsy))
        eddy_center_lat1,eddy_center_lon1 = utm.to_latlon(eddy_center_x[eddies-1],
                                        eddy_center_y[eddies-1],zone_eddy,zone_eddy_str)
        eddy_center_lat.append(eddy_center_lat1)
        eddy_center_lon.append(eddy_center_lon1)
        eddy_dir.append(winding_dir[k[0][0]])
        eddy_angular_vel.append(np.mean(eddy_angular_vel1)); # % avg angular velocity in deg/s
        eddy_streamlines.append(len(k));
       
        # %
        # %     [d,l] = ellipse_fit(ptsx,ptsy); % d is directions, l is lengths
        # % these three lines are the function
        M = np.cov(ptsx,ptsy);
        l,d = np.linalg.eig(M)
        dx1 = d[0][0]
        dy1 = d[1][0]
        l1 = np.sqrt(np.abs(l[0])); # % sqrt?
        l2 = np.sqrt(np.abs(l[1])); # % sqrt?
        
         
        ellipse_theta = np.arctan2(dy1,dx1)*180/np.pi
        eddy_ellipse_theta.append(ellipse_theta)
        eddy_length_x.append(l1); # % in meters
        eddy_length_y.append(l2);
           
        if debug_plot == 1: # % ellipse plots
            it = np.array(range(0,360))
            # it = 1:360;
            ellipse_x = l1*np.cos(it*np.pi/180);
            ellipse_y = l2*np.sin(it*np.pi/180);
            # % Create rotation matrix
            R = np.matrix([[np.cos(ellipse_theta*np.pi/180), -np.sin(ellipse_theta*np.pi/180)],
                          [np.sin(ellipse_theta*np.pi/180), np.cos(ellipse_theta*np.pi/180)]]);
            # % Rotate points
            eddypts = R*[ellipse_x , ellipse_y];
            ellipse_x = eddypts[0];
            ellipse_y = eddypts[1];
        
            ellipse_x1 = eddy_center_x[eddies-1] + ellipse_x
            ellipse_x1 = np.array(ellipse_x1.flatten())
            ellipse_y1 = eddy_center_y[eddies-1] + ellipse_y
            ellipse_y1 = np.array(ellipse_y1.flatten())
        
            eddy_center_lat_pts,eddy_center_lon_pts = utm.to_latlon(ellipse_x1,
                                        ellipse_y1,zone_eddy,zone_eddy_str)
            
            ptslat,ptslon = utm.to_latlon(ptsx,
                                        ptsy,zone_eddy,zone_eddy_str)
            
            fig = plt.figure(figsize=(10,6))
            plt.plot(eddy_center_lon1,eddy_center_lat1,'kx')#,'linewidth',4),'markersize',12)
            plt.plot(eddy_center_lon_pts,eddy_center_lat_pts,'b.')
            plt.plot(ptslon,ptslat,'r.')
            plt.savefig(fig_name + '_eddy_' + str(eddies) + '.png')
            plt.close(fig)
#            print(eddy_center_lon1)
#            print(eddy_center_lat1)
    
        
    # %% streamlines for each eddy we found
    if debug_info == 1:
        print(str(eddies) + ' eddies found')
        for i_eddy in range(0,eddies):
            print('eddy #' + str(i_eddy+1) + ' has ' 
                  + str(eddy_streamlines[i_eddy]) + ' streamlines')
    
    # show time of calculation
    elapsed = time.time() - timer1
    if debug_info == 1:
        print('processing time = ' + str(np.round(elapsed)) + ' seconds')   
    
    # convert to numpy arrays
    eddy_center_lat = np.array(eddy_center_lat) # eddy center lon/lat
    eddy_center_lon = np.array(eddy_center_lon)
    eddy_dir = np.array(eddy_dir) # eddy direction, math notation
    eddy_angular_vel = np.array(eddy_angular_vel) # angular velocty (deg/sec) math notation of eddy
    eddy_streamlines = np.array(eddy_streamlines) # number of streamlines for each eddy
    eddy_length_x = np.array(eddy_length_x) # ellipse fit x length (meters)
    eddy_length_y = np.array(eddy_length_y) # ellipse fit y length (meters)
    eddy_ellipse_theta = np.array(eddy_ellipse_theta) # rotation of eddy ellipse fit, math degrees
    # nc_time_ordinal = np.array(nc_time_ordinal)
    
    # open a netCDF file to write
    ncout = Dataset(ncout_name, 'w', format='NETCDF4')
    # define axis size, number of eddies
    ncout.createDimension('Number of eddies', eddies)  
    
    # create variable array
    eddy_center_lat_out = ncout.createVariable('eddy_center_lat', np.dtype('double').char, ('Number of eddies'))
    eddy_center_lat_out.long_name = ' eddy center lat'
    eddy_center_lat_out.units = 'Degrees'
    
    # create variable array
    eddy_center_lon_out = ncout.createVariable('eddy_center_lon', np.dtype('double').char, ('Number of eddies'))
    eddy_center_lon_out.long_name = ' eddy center lat'
    eddy_center_lon_out.units = 'Degrees'
    
    # create variable array
    eddy_dir_out = ncout.createVariable('eddy_dir', np.dtype('double').char, ('Number of eddies'))
    eddy_dir_out.long_name = 'Direction (math notation) of eddy rotation'
    eddy_dir_out.units = '+/-1'
    
    # create variable array
    eddy_angular_vel_out = ncout.createVariable('eddy_angular_vel', np.dtype('double').char, ('Number of eddies'))
    eddy_angular_vel_out.long_name = 'angular velocty (deg/sec) math notation of eddy'
    eddy_angular_vel_out.units = 'Degrees/sec'
    
    # create variable array
    eddy_streamlines_out = ncout.createVariable('eddy_streamlines', np.dtype('double').char, ('Number of eddies'))
    eddy_streamlines_out.long_name = 'number of streamlines for each eddy'
    eddy_streamlines_out.units = 'Integer'
    
    # create variable array
    eddy_length_x_out = ncout.createVariable('eddy_length_x', np.dtype('double').char, ('Number of eddies'))
    eddy_length_x_out.long_name = 'ellipse fit x length (meters)'
    eddy_length_x_out.units = 'meters'
    
    # create variable array
    eddy_length_y_out = ncout.createVariable('eddy_length_y', np.dtype('double').char, ('Number of eddies'))
    eddy_length_y_out.long_name = 'Dellipse fit y length (meters)'
    eddy_length_y_out.units = 'meters'
    
    # create variable array
    eddy_ellipse_theta_out = ncout.createVariable('eddy_ellipse_theta', np.dtype('double').char, ('Number of eddies'))
    eddy_ellipse_theta_out.long_name = 'rotation of eddy ellipse fit, math notation, degrees'
    eddy_ellipse_theta_out.units = 'Degrees'
    
    # # create variable array
    # nc_time_ordinal_out = ncout.createVariable('nc_time_ordinal', np.dtype('object').char, ('1'))
    # nc_time_ordinal_out.long_name = 'Python ordinal time in UTC'
    # nc_time_ordinal_out.units = 'Degrees'
    
    
    # copy axis from original dataset
    eddy_center_lat_out[:]      = eddy_center_lat[:]
    eddy_center_lon_out[:]      = eddy_center_lon[:]
    eddy_dir_out[:]             = eddy_dir[:]
    eddy_angular_vel_out[:]     = eddy_angular_vel[:]
    eddy_streamlines_out[:]     = eddy_streamlines[:]
    eddy_length_x_out[:]        = eddy_length_x[:]
    eddy_length_y_out[:]        = eddy_length_y[:]
    eddy_ellipse_theta_out[:]   = eddy_ellipse_theta[:]
    # nc_time_ordinal_out         = nc_time_ordinal[:]
    
    ncout.close()
    
    if debug_info == 1:
        print('analysis on this file completed')
    
    # if plt_debug == 1 % save plot
    #     savefig(sn)
    #     print(sn,'-dpng')
    #     if numfiles ~= 1
    #         close(h)
    #     end
    # end
        
    # break
    
    
    # %% 
    # break
    if file_count == numfiles or file_count > numfiles:
        break # only one file for now  
        
        
        




