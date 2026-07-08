%% lon lat limits
params.lonlatlimits = 0; % (=1) enable limiting to subregion, set 0 to do entire region

%% parameters

params.plt_debug            = 1;    % debug plots
params.plt_debug_streams    = 0;    % streamline plots
params.save_streams         = 1;    % save lat lon of streamlines
params.plt_quiver           = 1;    % plot quiver
params.plt_debug_sv         = 1;    % save debug plot
params.highrescoast         = 1;    % high resolution coastline
params.bathy                = 1;    % plot bathy levels
params.levels = [-10 -25 -50 -100 -250 -500 -1000 -1500 -2000 -2500 -3000 -4000 -5000]; % isobaths 
params.plt_close            = 1;    % close figures after each timestep
params.print_res            = 300;  % resolution of plot
params.plt_each_eddy        = 1;    % zooms and plots figure for each eddy
params.plt_min_streams      = 5;

params.winding_thres       = 300;   % how much winding (in degrees) to be a eddy
params.d_thres             = 10;    % streamline starting and ending max distance in km
params.baddir_thres        = 15;    % breaks at 15 degrees of bad dirs
% params.param_center_dist   = 5;     % distance between streamlines centers
params.new_dist_thres      = 1;     % Cahl et al., 2023 method of clustering (no param_center_dist_used)

params.param_center_dist   = 2*sqrt(2)*radar_km_resolution; % distance between streamlines centers
params.radar_km_resolution = radar_km_resolution;

params.dij = 100; % streamlines in areas of lons(i-dij:i+dij,j-dij:j+dij)
% for 1km resolution, dij = 100 for max diameter streamline of 1-200 km