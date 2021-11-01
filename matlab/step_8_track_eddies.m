function step_8_track_eddies(bool_avg)

%%
cd 'C:\Users\doug\Google Drive\MATLAB\current work'

if bool_avg == 1
    %     infile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_25_hr_avg_eddies_data.mat';
%     outfile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_25_hr_avg_eddies_data_tracks.mat';
    infile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_25_hr_avg_eddies_data_no_streams.mat';
    outfile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_25_hr_avg_eddies_data_no_streams_tracks.mat';
end

if bool_avg == 0
    %     infile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_eddies_data.mat';
%     outfile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_eddies_data_tracks.mat';
    infile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_eddies_data_no_streams.mat';
    outfile = 'LNG_joined_non_QAQC_mask_1_2012_5_2014_7_wildpt_removed_eddies_data_no_streams_tracks.mat';
    
end

a = load(infile);
%%

% tracks
dist_threshold = 4000; % 4 km
skip_threshold = 5; % number of timesteps away where we continue the eddy track. 2 allows 1 missing timestep, i.e. look two timesteps back


%%
j = 0;
for i = 1:length(a.direction)
    j = j + length(a.direction{i});
    jj(i) = length(a.direction{i});
end

disp(['Total eddies = ' num2str(j)])
disp(['Max eddies in a single timestep = ' num2str(max(jj))])
%%


it = 1; % timestep 1
i_cluster = length(a.x_cluster_center{it}); % i_eddy is the cluster

for i = 1:i_cluster
    x_center{it} = a.x_cluster_center{it}(i);
    y_center{it} = a.y_cluster_center{it}(i);
    vec1x{it} = a.vec1x{it}(i);
    vec1y{it} = a.vec1y{it}(i);
    vec2x{it} = a.vec2x{it}(i);
    vec2y{it} = a.vec2y{it}(i);
    eig1{it} = a.eig1{it}(i);
    eig2{it} = a.eig2{it}(i);
    omega{it} = a.omega{it}(i); % rads per sec
    direction{it} = a.direction{it}(i); % +1 is cyclonic, -1 is anticyclonic
    num_streams{it} = a.num_streams{it}(i); % intensity
    Time{it} = a.Time(it);
    timegap(it) = 1;
end
%%

total_tracks = i_cluster;
disp('Starting track finding routine')
disp(['Total number of timestep to analyze = ' num2str(length(a.Time))])
pause(2)

for it = 2:length(a.Time);
    disp(['Timestep = ' num2str(it)])
    i_cluster = length(a.x_cluster_center{it}); % number of eddies at this timestep
    %     i_cluster = num_clusters(it);
    for i=1:i_cluster
        curr_x_center = a.x_cluster_center{it}(i);
        curr_y_center = a.y_cluster_center{it}(i);
        
        %         currentpos = [x(it,i) y(it,i)];
        %         cur_size = [eigvals1(it,i) eigvals2(it,i)];
        
        % check over all tracks with timesteps from before
        succ = 0;
        i1_start = find(timegap <= skip_threshold,1);
        for i1=i1_start:total_tracks;
          
            dist = sqrt( (x_center{i1}(end) - curr_x_center)^2 + (y_center{i1}(end) - curr_y_center)^2 );
            %             dist = norm([tracks(i1).center(end,1) tracks(i1).center(end,2)]  - currentpos);
            if dist < dist_threshold
                for i_tgap = 1:skip_threshold
                    if timegap(i1) <= i_tgap
                        x_center{i1}    = [x_center{i1} ; a.x_cluster_center{it}(i)];
                        y_center{i1}    = [y_center{i1} ; a.y_cluster_center{it}(i)];
                        vec1x{i1}       = [vec1x{i1}; a.vec1x{it}(i)];
                        vec1y{i1}       = [vec1y{i1} ; a.vec1y{it}(i)];
                        vec2x{i1}       = [vec2x{i1}; a.vec2x{it}(i)];
                        vec2y{i1}       = [vec2y{i1} ; a.vec2y{it}(i)];
                        eig1{i1}        = [eig1{i1} ; a.eig1{it}(i)];
                        eig2{i1}        = [eig2{i1} ; a.eig2{it}(i)];
                        omega{i1}       = [omega{i1} ; a.omega{it}(i)];
                        direction{i1}   = [direction{i1}; a.direction{it}(i)]; % +1 is cyclonic, -1 is anticyclonic
                        num_streams{i1} = [num_streams{i1} ; a.num_streams{it}(i)]; % intensity
                        Time{i1}        = [Time{i1} ; a.Time(it)];
%                         disp(['Time gaps away = ' num2str(timegap(i1))])
                        timegap(i1)     = 0; % zero because it's gonna add one right after this loop
                        succ = 1;
                        break;
                    end
                end
%                 if tracks(i1).timegap == 1; % previous timestep
%                     % disp('continuation')
%                     tracks(i1).center = [tracks(i1).center ; currentpos];
%                     tracks(i1).time = [tracks(i1).time ; it];
%                     tracks(i1).length = tracks(i1).length + 1;
%                     tracks(i1).timegap = 0;
%                     tracks(i1).size = [tracks(i1).size ; cur_size];
%                     tracks(i1).omega = [tracks(i1).omega ; w(it,i)];
%                     succ = 1;
%                     break;
%                 elseif tracks(i1).timegap == 2; % 2 timesteps before
%                     % disp('missed just 1 frame')
%                     tracks(i1).center = [tracks(i1).center ; currentpos];
%                     tracks(i1).time = [tracks(i1).time ; it];
%                     tracks(i1).length = tracks(i1).length + 1;
%                     tracks(i1).timegap = 0;
%                     tracks(i1).size = [tracks(i1).size ; cur_size];
%                     tracks(i1).omega = [tracks(i1).omega ; w(it,i)];
%                     succ = 1;
%                     break;
%                 end
            end
        end
        % if not close to others, create new track
        if succ == 0
            total_tracks = total_tracks + 1;
            x_center{total_tracks} = a.x_cluster_center{it}(i);
            y_center{total_tracks} = a.y_cluster_center{it}(i);
            vec1x{total_tracks} = a.vec1x{it}(i);
            vec1y{total_tracks} = a.vec1y{it}(i);
            vec2x{total_tracks} = a.vec2x{it}(i);
            vec2y{total_tracks} = a.vec2y{it}(i);
            eig1{total_tracks} = a.eig1{it}(i);
            eig2{total_tracks} = a.eig2{it}(i);
            omega{total_tracks} = a.omega{it}(i);
            direction{total_tracks} = a.direction{it}(i); % +1 is cyclonic, -1 is anticyclonic
            num_streams{total_tracks} = a.num_streams{it}(i); % intensity
            Time{total_tracks} = a.Time(it);
            timegap(total_tracks) = 0; % zero because it's gonna add one right after this loop
%             tracks(total_tracks).center = [x(it,i) y(it,i)];
%             tracks(total_tracks).time = it;
%             tracks(total_tracks).length = 1;
%             tracks(total_tracks).timegap = 0;
%             tracks(total_tracks).size = [eigvals1(it,i) eigvals2(it,i)];
%             tracks(total_tracks).year = str2num(year);
%             tracks(total_tracks).month = str2num(month);
%             tracks(total_tracks).omega = w(it,i);
        end
    end
    
    % add 1 timestep to the counter
    for i1=i1_start:total_tracks;
        timegap(i1) = timegap(i1) + 1;
%         tracks(i1).timegap = tracks(i1).timegap + 1;
    end
    
end
%%
% save(['C:\Users\doug\Documents\MATLAB\eddy\eddy_data_' year '_' month '.mat'],'tracks','total_tracks','dist_threshold','-append')
% disp(['saved ' 'C:\Users\doug\Documents\MATLAB\eddy\eddy_data_' year '_' month '.mat'])

% save(outfile,'tracks','total_tracks','dist_threshold','skip_threshold')
save(outfile,'total_tracks','x_center','y_center','vec1x','vec1y','vec2x','vec2y','eig1','eig2','omega','direction','num_streams','Time','timegap','dist_threshold','skip_threshold')


end





























