%% to run your own data set (lon, lat, u, v, time)
%
% Douglas Cahl, PhD 
% University of South Carolina 2023


clear
%% load data
load data/data2.mat
% u1(:,:101) v1(:,:101) time1(101) lons(:,:) lats(:,:)
        
%% lon lat limits
params.lonlatlimits = 0; % set =1 to enable limiting to subregion, set 0 to do entire region

%% parameters
% params.save_name = fn;

params.plt_debug            = 0;    % debug plots
params.plt_debug_streams    = 0;    % streamline plots
params.plt_debug_sv         = 0;    % save debug plot
params.highrescoast         = 0;    % high resolution coastline
params.bathy                = 0;    % plot bathy levels
params.levels = [-10 -25 -50 -100 -250 -500 -1000 -1500 -2000 -2500 -3000 -4000 -5000]; % isobaths 
params.plt_close            = 0;    % close figures after each timestep
params.print_res            = 300;  % resolution of plot

params.winding_thres       = 300;   % how much winding (in degrees) to be a eddy
params.d_thres             = 10;    % streamline starting and ending max distance in km
params.baddir_thres        = 15;    % breaks at 15 degrees of bad dirs
params.param_center_dist   = 5;     % only used if new_dist_thres=0; distance between streamlines centers
params.new_dist_thres      = 1;     % Cahl et al., 2023 method of clustering (no param_center_dist_used)



%% start running data
numfiles = -1; % number of files to run

counter = 0;
for i = 1:length(time1)
    disp(i)
    counter = counter + 1;

    % data
    time = time1(i);
    u = u1(:,:,i);
    v = v1(:,:,i);
    
    % results filenames
    params.save_name = ['data/results/data2_' num2str(i)];

    % call eddy finding
    eddy_subroutine(params,u,v,time,lons,lats)
    
    % run numfiles
    if counter >= numfiles && numfiles > 0
        break
    end


end





