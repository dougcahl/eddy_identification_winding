%% to run your own data set (lon, lat, u, v, time)
%
% Douglas Cahl, PhD 
% University of South Carolina 2023

%% parameters
params_eddy1;

%% start processing
file_count = 0;
while true
    nc_time = nc_time_start + file_count; % hours since 2019-09-29 20:00:00.000 UTC
    if nc_time > tnum
        break
    end

    timer1a = tic;
    file_count = file_count + 1;
   
    disp(['running analysis on file ' num2str(nc_time)])
   
    % url name stuff, ex. fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'
    if avg_25hr == 1
        uname = 'u_mean';
        vname = 'v_mean';
        if (region_num == 6) && (radar_km_resolution == 6)
            uname = 'u';
            vname = 'v';
        end
    else
        uname = 'u';
        vname = 'v';
    end

    % opendap data
    nc = ['http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/' region ...
    '/' num2str(radar_km_resolution) ...
    'km/' hr_name ...
    '/RTV/' fullregionname ...
    '_' num2str(radar_km_resolution) ...
     'km_Resolution_' hr_name_full ...
    '_RTV_best.ncd?' 'lat[0:1:' num2str(latnum-1) ...
    '],lon[0:1:' num2str(lonnum-1) ...
    '],time[' num2str(nc_time) ...
    '],' uname ...
    '[' num2str(nc_time) ...
    '][0:1:' num2str(latnum-1) ...
    '][0:1:' num2str(lonnum-1) ...
    '],' vname ...
    '[' num2str(nc_time) ...
    '][0:1:' num2str(latnum-1) ...
    '][0:1:' num2str(lonnum-1) ...
    ']'];

%     disp(url)
%     url = 'http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USEGC/6km/25hr/RTV/HFRADAR_US_East_and_Gulf_Coast_6km_Resolution_25_Hour_Average_RTV_best.ncd?lat[0:1:459],lon[0:1:700],time[11746],u_mean[11746][0:1:459][0:1:700],v_mean[11746][0:1:459][0:1:700]'
    % region = 'USEGC'
    % nc_time = 11764 % hours since 2019-09-29 20:00:00.000 UTC
    
    try
        disp('downloading data ...')
        disp(nc)
        time = ncread(nc,'time');
        lon = double(ncread(nc,'lon'));
        lat = double(ncread(nc,'lat'));
        if avg_25hr == 1
            if (region_num == 6) && (radar_km_resolution == 6)
                u = ncread(nc,'u');
                v = ncread(nc,'v');
            else
                u = ncread(nc,'u_mean');
                v = ncread(nc,'v_mean');
            end
        else
            u = ncread(nc,'u');
            v = ncread(nc,'v');
    
            u = ncread(nc,'u');
            v = ncread(nc,'v');
        end
        
        disp('downloaded')
    catch
        disp('download failed, file does not exist or website down')
        disp(nc)
        break
    end
    
    % mesh lon/lat
    [lats,lons] = meshgrid(lat,lon);
    

    % results filenames
    params.save_data_name = [dout_data region '_' num2str(radar_km_resolution) 'km_' num2str(nc_time)];
    params.save_fig_name  = [dout_fig region '_' num2str(radar_km_resolution) 'km_' num2str(nc_time)];

    % call eddy finding
    disp('finding eddies...')
    timer1a = tic;
    eddy_subroutine(params,u,v,time,lons,lats)
    disp(['analysis complete:' num2str(toc(timer1a)) ' seconds'])

  
        
    % run numfiles
    if file_count == numfiles || file_count > numfiles
        break
    end
end








