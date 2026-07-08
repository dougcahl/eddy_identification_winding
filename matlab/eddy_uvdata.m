%% to run your own data set (lon, lat, u, v, time)
%
% Douglas Cahl, PhD 
% University of South Carolina 2023


clear
%% load data
data_dir = 'data/';
fn = '202004030900_hfr_usegc_2km_rtv_uwls_NDBC'; % filename
ffn = [data_dir fn];
nc = [ffn '.nc']; % load from netcdf
lon = double(ncread(nc,'lon'));
lat = double(ncread(nc,'lat'));
u = ncread(nc,'u');
v = ncread(nc,'v');
time = ncread(nc,'time');
[lats,lons] = meshgrid(lat,lon);
        
%% lon lat limits
params.lonlatlimits = 1; % (=1) enable limiting to subregion, set 0 to do entire region
% lat_min = 38.5;
% lat_max = 39;
% lon_min = -75.3;
% lon_max = -74.6;

params.lat_min = 38.5;
params.lat_max = 39.4;
params.lon_min = -75.5;
params.lon_max = -74.5;



%% parameters
params.save_data_name = [data_dir 'results/' fn];
params.save_fig_name = [data_dir 'results_figs/' fn];

params.plt_debug            = 1;    % debug plots
params.plt_debug_streams    = 0;    % streamline plots
params.save_streams         = 1;    % save lat lon of streamlines
params.plt_quiver           = 1;    % plot quiver
params.plt_debug_sv         = 1;    % save debug plot
params.highrescoast         = 1;    % high resolution coastline
params.bathy                = 1;    % plot bathy levels
params.levels = [-10 -25 -50 -100 -250 -500 -1000 -1500 -2000 -2500 -3000 -4000 -5000]; % isobaths 
params.plt_close            = 0;    % close figures after each timestep
params.print_res            = 300;  % resolution of plot
params.plt_each_eddy        = 0;    % zooms and plots figure for each eddy
params.plt_min_streams      = 5;

params.winding_thres       = 300;   % how much winding (in degrees) to be a eddy
params.d_thres             = 10;    % streamline starting and ending max distance in km
params.baddir_thres        = 15;    % breaks at 15 degrees of bad dirs
radar_km_resolution        = 2;
params.param_center_dist   = 2*sqrt(2)*radar_km_resolution;  % distance between streamlines centers
params.new_dist_thres      = 1;     % Cahl et al., 2023 method of clustering (no param_center_dist_used)
params.radar_km_resolution = radar_km_resolution;

params.dij = 50; % streamlines in areas of lons(i-dij:i+dij,j-dij:j+dij)
% for 1km resolution, dij = 100 for max diameter streamline of 1-200 km

%% call eddy finding
eddy_subroutine(params,u,v,time,lons,lats,lon,lat)





