%% master for eddy identification, tracking and some eddy track analysis
%
% Douglas Cahl, PhD 
% University of South Carolina 2023


%% single eddy identificaiton from netcdf file
eddy_uvdata;

%% identify eddies in consecutive time steps
% runs all data in data/data2.mat, change params to suppress debug plots
eddy_uvdata_loop;

%% eddy tracking
% tracks the eddies found from above script
eddy_tracking;

%% plots the eddy tracks found from above script and some data analysis
analyze_eddy_tracks;