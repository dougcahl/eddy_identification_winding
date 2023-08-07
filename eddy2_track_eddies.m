%% eddy tracking 
% if not on first timestep
% check each eddy on this timestep for distance from previous eddies
% if distance < eddy_track_dist_param
% andif same winding dir as before
% then this eddy is a continuation of the last eddy
% check if this eddy track exists
% if it does, add new data point
% else look back eddy_track_time_param hours/timesteps for eddy
% else if no previous eddy, make new eddy track
%
% Douglas Cahl, PhD 
% University of South Carolina 2023


%% parameters
params_eddy2;


%% directory setup
name_dir = dout_data; % results folder
name_pre  = [region '_' num2str(radar_km_resolution) 'km_']; % name prefix, e.g., data2_1.mat , data2_2.mat , ...
outfile = [tracks_dir name_pre 'tracks'];



%% load data
j = 0;
f = dir(name_dir);
for i = 3:length(f)
    fn = f(i).name;
    if ~strcmp(fn(end-3:end),'.mat') || ~contains(fn,name_pre) || contains(fn,'tracks')
        continue
    end
    disp(fn)
    data = load([name_dir fn]);
    if ~isempty(data.eddy_streamlines_lat)
        j = j + 1;
        tmp.eddy_angular_vel{j}   = data.eddy_angular_vel;
        tmp.eddy_center_lat{j}    = data.eddy_center_lat;
        tmp.eddy_center_lon{j}    = data.eddy_center_lon;
        tmp.eddy_dir{j}           = data.eddy_dir;
        tmp.eddy_ellipse_theta{j} = data.eddy_ellipse_theta;
        tmp.eddy_length_x{j}      = data.eddy_length_x;
        tmp.eddy_length_y{j}      = data.eddy_length_y;
        tmp.eddy_streamlines{j}   = data.eddy_streamlines;
    
        tmp.time(j) = data.time;
    end
    
    clear data
end

if ~exist('tmp','var')
    return
end

%% sort by time
[time,timei] = sort(tmp.time);
for i = 1:length(timei)
    j = timei(i);
    eddy_angular_vel{i} = tmp.eddy_angular_vel{j};
    eddy_center_lat{i} = tmp.eddy_center_lat{j};
    eddy_center_lon{i} =  tmp.eddy_center_lon{j};
    eddy_dir{i} =  tmp.eddy_dir{j};
    eddy_ellipse_theta{i} =  tmp.eddy_ellipse_theta{j};
    eddy_length_x{i} =  tmp.eddy_length_x{j};
    eddy_length_y{i} =  tmp.eddy_length_y{j};
    eddy_streamlines{i} =  tmp.eddy_streamlines{j};
end
clear tmp

%%
j = 0;
for i = 1:length(eddy_dir)
    j = j + length(eddy_dir{i});
    jj(i) = length(eddy_dir{i});
    mxstreams = max(eddy_streamlines{i});
    if ~isempty(mxstreams)
        jjstreams(i) = max(eddy_streamlines{i});
    else
        jjstreams(i) = 0;
    end
end

disp(['Total eddies = ' num2str(j)])
disp(['Max eddies in a single timestep = ' num2str(max(jj))])
disp(['Max streamlines for an eddy = ' num2str(max(jjstreams))])


%% no tracking for first timestep
it = 1; % timestep 1
i_cluster = length(eddy_center_lat{it}); % i_eddy is the cluster
for i = 1:i_cluster
    lat_center{i} = eddy_center_lat{it}(i);
    lon_center{i} = eddy_center_lon{it}(i);
    ellipse_theta{i} = eddy_ellipse_theta{it}(i);
    eig1{i} = eddy_length_x{it}(i);
    eig2{i} = eddy_length_y{it}(i);
    omega{i} = eddy_angular_vel{it}(i); % rads per sec
    direction{i} = eddy_dir{it}(i); % +1 is cyclonic, -1 is anticyclonic
    num_streams{i} = eddy_streamlines{it}(i); % intensity
    Time{i} = time(it);
    timegap(i) = 1;
end
if i_cluster == 0
    timegap = 1;
end

%%
i1_start = 1;
total_tracks = i_cluster;
disp('Starting track finding routine')
disp(['Total number of timestep to analyze = ' num2str(length(time))])
% disp(['Timestep = ' num2str(it)])

for it = 2:length(time)
%     disp(['Timestep = ' num2str(it)])
    i_cluster = length(eddy_center_lat{it}); % number of eddies at this timestep
    for i=1:i_cluster
        curr_lat_center = eddy_center_lat{it}(i);
        if isnan(curr_lat_center)
            continue
        end
        curr_lon_center = eddy_center_lon{it}(i);
        curr_ellipse_x = eddy_length_x{it}(i);
        curr_ellipse_y = eddy_length_y{it}(i);
        curr_ellipse_th = eddy_ellipse_theta{it}(i);

        % check over all tracks with timesteps from before
        succ = 0;
        i1_start = find(timegap <= eddy_track_time_param,1);
        for i1=i1_start:total_tracks
            
            eddy_continuation = 0; % 
            if new_dist_thres % Cahl et al., 2023 method (if center is within another streamline, they are the same eddy)
                % in km
                i_ellipse = 30:30:360;
                ellipse_x = eig1{i1}(end)*cosd(i_ellipse);
                ellipse_y = eig2{i1}(end)*sind(i_ellipse);
                ellipse_th = ellipse_theta{i1}(end);
                % Create rotation matrix
                R = [cosd(ellipse_th) -sind(ellipse_th); sind(ellipse_th) cosd(ellipse_th)];
                % Rotate points
                eddypts = R*[ellipse_x ; ellipse_y];
                ellipse_x = eddypts(1,:);
                ellipse_y = eddypts(2,:);
                
                 % in km
                 [x,y,UTMzone] = geog2utm_nodisp(lon_center{i1}(end),lat_center{i1}(end),curr_lon_center,curr_lat_center);
                ellipse_x = ellipse_x + x;
                ellipse_y = ellipse_y + y;

                % is old center within new eddy ellipse
                in1 = inpolygon(0,0,ellipse_x,ellipse_y);


                % in km
                i_ellipse = 30:30:360;
                % Create rotation matrix
                R = [cosd(curr_ellipse_th) -sind(curr_ellipse_th); sind(curr_ellipse_th) cosd(curr_ellipse_th)];
                % Rotate points
                eddypts = R*[curr_ellipse_x ; curr_ellipse_y];
                ellipse_x1 = eddypts(1,:);
                ellipse_y1 = eddypts(2,:);

                % is new center within old eddy ellipse
                in2 = inpolygon(x,y,ellipse_x1,ellipse_y1);


                if in1 || in2 % 
                    eddy_continuation = 1;
                else
                    % in km
                    [x,y,UTMzone] = geog2utm_nodisp(lon_center{i1}(end),lat_center{i1}(end),curr_lon_center,curr_lat_center);
                    dist = sqrt(x^2+y^2);
                    % if centers are close and winding direction is the same
                    if dist < eddy_track_dist_param
                        eddy_continuation = 1;
                    end
                end
        
            else % Sadarjoen and Post 2000 method

                % in km
                [x,y,UTMzone] = geog2utm_nodisp(lon_center{i1}(end),lat_center{i1}(end),curr_lon_center,curr_lat_center);
                dist = sqrt(x^2+y^2);
                if dist < eddy_track_dist_param
                    eddy_continuation = 1;
                    
                end
            end
            if eddy_continuation
                for i_tgap = 1:eddy_track_time_param
                    if timegap(i1) <= i_tgap
                        lat_center{i1}    = [lat_center{i1} ; eddy_center_lat{it}(i)];
                        lon_center{i1}    = [lon_center{i1} ; eddy_center_lon{it}(i)];
                        ellipse_theta{i1}       = [ellipse_theta{i1} ; eddy_ellipse_theta{it}(i)];
                        eig1{i1}        = [eig1{i1} ; eddy_length_x{it}(i)];
                        eig2{i1}        = [eig2{i1} ; eddy_length_y{it}(i)];
                        omega{i1}       = [omega{i1} ; eddy_angular_vel{it}(i)];
                        direction{i1}   = [direction{i1}; eddy_dir{it}(i)]; % +1 is cyclonic, -1 is anticyclonic
                        num_streams{i1} = [num_streams{i1} ; eddy_streamlines{it}(i)]; % intensity
                        Time{i1}        = [Time{i1} ; time(it)];
                        timegap(i1)     = 0; % zero because it's gonna add one right after this loop
                        succ = 1;
                        break;
                    end
                end
            end

        end
        % if not close to others, create new track
        if succ == 0
            total_tracks = total_tracks + 1;
            lat_center{total_tracks} = eddy_center_lat{it}(i);
            lon_center{total_tracks} = eddy_center_lon{it}(i);
            ellipse_theta{total_tracks} = eddy_ellipse_theta{it}(i);
            eig1{total_tracks} = eddy_length_x{it}(i);
            eig2{total_tracks} = eddy_length_y{it}(i);
            omega{total_tracks} = eddy_angular_vel{it}(i);
            direction{total_tracks} = eddy_dir{it}(i); % +1 is cyclonic, -1 is anticyclonic
            num_streams{total_tracks} = eddy_streamlines{it}(i); % intensity
            Time{total_tracks} = time(it);
            timegap(total_tracks) = 0; % zero because it's gonna add one right after this loop
        end
    end
    % add 1 timestep to the counter
    for i1=i1_start:total_tracks
        timegap(i1) = timegap(i1) + 1;
    end
end
disp('Finished eddy track finding')

%% save here
save(outfile,'total_tracks','lat_center','lon_center','ellipse_theta',...
    'eig1','eig2','omega','direction','num_streams','Time','timegap',...
    'eddy_track_dist_param','eddy_track_time_param')


