%% parameters
eddy_track_dist_param = sqrt(2)*radar_km_resolution; % distance in km from previous eddy center for eddy to be continuation
eddy_track_time_param = 4;  % how many timesteps can there be between identifications
new_dist_thres        = 1;

