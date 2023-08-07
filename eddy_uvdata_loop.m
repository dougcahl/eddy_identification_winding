%% to run your own data set (lon, lat, u, v, time)
%
% Douglas Cahl, PhD 
% University of South Carolina 2023


clear
%% load data
load data/data2.mat
% u1(:,:101) v1(:,:101) time1(101) lons(:,:) lats(:,:)
%% if from WERA site need to flip so that Lats increase with increasing index
lats = fliplr(lats);
lons = fliplr(lons);
u1 = fliplr(u1);
v1 = fliplr(v1);
%% lon lat limits
params.lonlatlimits = 0; % set =1 to enable limiting to subregion, set 0 to do entire region

%% parameters
% params.save_name = fn;

params.plt_debug            = 0;    % debug plots
params.plt_debug_streams    = 0;    % streamline plots, doesn't work with plt_each_eddy=1
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
radar_km_resolution        = 3;
params.param_center_dist   = 2*sqrt(2)*radar_km_resolution;  % distance between streamlines centers
params.new_dist_thres      = 1;     % Cahl et al., 2023 method of clustering
params.radar_km_resolution = radar_km_resolution;

params.dij = 50; % streamlines in areas of lons(i-dij:i+dij,j-dij:j+dij)
% for 1km resolution, dij = 100 for max diameter streamline of 1-200 km

%% start running data
numfiles = -1; % number of files to run

counter = 0;
% for i = 61 % has three eddies
for i = 1:length(time1) 
    disp(i)
    counter = counter + 1;

    % data
    time = time1(i);
    u = u1(:,:,i);
    v = v1(:,:,i);
    
    % results filenames
    params.save_data_name = ['data/results/data2_' num2str(i)];
    params.save_fig_name = ['data/results_figs/data2_' num2str(i)];

    % call eddy finding
    disp('finding eddies...')
    timer1a = tic;
    eddy_subroutine(params,u,v,time,lons,lats)
    disp(['analysis complete:' num2str(toc(timer1a)) ' seconds'])
    
    % run numfiles
    if counter >= numfiles && numfiles > 0
        break
    end


end





