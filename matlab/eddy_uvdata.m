%% to run your own data set (lon, lat, u, v)
clear
%% load uv data, this is just the data from the netcdf in the next section
% load data1.mat

%% or load from netcdf
nc = '202004030900_hfr_usegc_2km_rtv_uwls_NDBC.nc'; % this one has some eddies in it
lon = double(ncread(nc,'lon'));
lat = double(ncread(nc,'lat'));
u = ncread(nc,'u');
v = ncread(nc,'v');
time = ncread(nc,'time');
        
%% lon lat limits
lonlatlimits = 1; % (=1) enable limiting to subregion, set 0 to do entire region
lat_min = 38.5;
lat_max = 39;
lon_min = -75.3;
lon_max = -74.6;
highrescoast = 0; % high resolution coastline for smaller regions

%% parameters
plt_debug = 1; % debug plots
plt_debug_streams = 0; % streamline plots

numfiles = 1; % if not =1 it will close figures after running

radar_km_resolution = 2;    % 1, 2 or 6 km nc file
res_x_mult          = 1;    % for higher resolution, if >1 results in this squared processing time (i.e. 3 => 9 times longer)
winding_thres       = 270;  % how much winding (in degrees) to be a eddy
d_thres             = .25;   % streamline starting and ending max distance in km
baddir_thres        = 15;   % breaks at 15 degrees of bad dirs
param_center_dist   = .25;   % how many km are clusters together

%% output filenames
% figure ouput filename
sn = 'D:\secoora_data_challenge\project\data_challenge_2020\eddy_analysis\202004010000_hfr_usegc_2km_rtv_uwls_NDBC';
% .mat output filename
snmat = 'D:\secoora_data_challenge\project\data_challenge_2020\eddy_analysis\mat_2km\202004010000_hfr_usegc_2km_rtv_uwls_NDBC';


%% call eddy finding
eddy_subroutine