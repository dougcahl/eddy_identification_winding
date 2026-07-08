%% parameters
min_steps = 6;      % minimum track length to plot (set 1 for all eddies)
radar_dt = 0.5;     % timestep in hours (used for histogram plot)

plt_individual = 1; % separate plot for each eddy track
% dlon = 1;           % longitude extent of plot is 2*dlon    
% dlat = 1;           % latitude extent of plot is 2*dlon
dlon = radar_km_resolution/6;
dlat = dlon;

highrescoast = 1;   % high resolution coastlines
print_res = 300;    % dpi of saved figure
plt = 1;            % plot and save figures
plt_close = 1;      % close after plotting 
lw = 1;             % linewidth
