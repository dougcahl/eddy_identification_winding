%% eddy finding script
clear
%% lon lat limits
lonlatlimits = 1; % (=1) enable limiting to subregion, set 0 to do entire region
% lat_min = 30;
% lat_max = 45;
% lon_min = -90;
% lon_max = -60;

lat_min = 38.5;
lat_max = 39;
lon_min = -75.3;
lon_max = -74.6;
highrescoast = 0; % high resolution coastline for smaller regions

%%
plt_debug = 1; % debug plots
plt_debug_streams = 0; % streamline plots

numfiles = 1; % how many files to run. Set to 0 to run all files

radar_km_resolution = 2;    % 1, 2 or 6 km nc file
res_x_mult          = 1;    % for higher resolution, if >1 results in this squared processing time (i.e. 3 => 9 times longer)
winding_thres       = 270;  % how much winding (in degrees) to be a eddy
d_thres             = .25;   % streamline starting and ending max distance in km
baddir_thres        = 15;   % breaks at 15 degrees of bad dirs
param_center_dist   = .25;   % how many km are clusters together

% if you have lots of data
data_dir = 'D:\secoora_data_challenge\project\data_challenge_2020\data\';
save_dir = 'D:\secoora_data_challenge\project\data_challenge_2020\eddy_analysis\';

numfiles_i = 1;
f = dir(data_dir);
for i = 3:length(f)
    fn = f(i).name; % netcdf filename
    fn_res = str2double(fn(24)); % 1, 2 or 6 km
    if fn_res == radar_km_resolution
        tic
        disp(fn)
        nc = [data_dir fn]; % netcdf full filename
%         nc = '202004030900_hfr_usegc_2km_rtv_uwls_NDBC.nc'; % this one has some eddies in it
        % save filename
        sn = [save_dir fn(1:end-3)];
        snmat = [save_dir 'mat_' num2str(radar_km_resolution) 'km\' fn(1:end-3)];
        
        lon = double(ncread(nc,'lon'));
        lat = double(ncread(nc,'lat'));
        u = ncread(nc,'u');
        v = ncread(nc,'v');
        time = ncread(nc,'time');
        
        
        % function params, lon, lat, u, v, sn, snmat
        
        
        %%
        [lats,lons] = meshgrid(lat,lon);
        
        if lonlatlimits == 1
            j = lons >= lon_min & lons <= lon_max & lats >= lat_min & lats <= lat_max; 
            u = u(j);
            v = v(j);
            lons = lons(j);
            lats = lats(j);
            
            j = lon >= lon_min & lon <= lon_max;
            lon1 = lon(j);
            j = lat >= lat_min & lat <= lat_max; 
            lat1 = lat(j);
            [lats1,lons1] = meshgrid(lat1,lon1);
            lons = reshape(lons,size(lats1));
            lats = reshape(lats,size(lats1));
            u = reshape(u,size(lats1));
            v = reshape(v,size(lats1));
            
        else
            lon_min = double(floor(min(lon)));
            lon_max = double(ceil(max(lon)));
            lat_min = double(floor(min(lat)));
            lat_max = double(ceil(max(lat)));
        end
        
        % higher res, this is updated to doing this in x,y km space in python script
        if res_x_mult > 1
            lon_n = length(lon);
            lat_n = length(lat);
            lon_high = linspace(lon_min,lon_max,lon_n*res_x_mult);
            lat_high = linspace(lat_min,lat_max,lat_n*res_x_mult);
            [lats_high,lons_high] = meshgrid(lat_high,lon_high);
        end

        %%  
        if plt_debug == 1
%             figure
%             m_quiver(lons,lats,u,v)
%             axis equal
           
            h = figure('color','w','renderer','painters');
            m_proj('Mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max])
            % m_gshhs_h('patch',[0.8 0.8 0.8]); % coastline
            if lonlatlimits == 1 && highrescoast == 1
                m_gshhs_i('patch',[0.8 0.8 0.8])
            else
                m_gshhs_l('patch',[0.8 0.8 0.8])
            end
                
            % levels = [-10:-5:-30 -40 -50 -100 -150 -200]; % isobaths
            % levels = [-100 -1000]; % isobaths
            % m_etopo2('contour',levels,'edgecolor','k','linewidth',1,'linestyle','--','ShowText','on');
            m_grid(); %'fontsize',14); % nice grid
            hold on
            m_quiver(lons,lats,u,v)
        end
        
        i_streams = 0;
        % convert lon lat to utm x,y in km
        ideg = 1;
        lon_deg = ceil(lon_max - lon_min);
        lat_deg = ceil(lat_max - lat_min);
        degmax = max(lon_deg,lat_deg);
        for ideg_lon = 0:lon_deg
            for ideg_lat = 0:lat_deg
                ill = (lons >= lon_min + ideg_lon) & ...
                    (lons <= lon_min + ideg_lon + 1) & ... % 1 deg increment
                    (lats >= lat_min + ideg_lat) & ...
                    (lats <= lat_min + ideg_lat + 1);
                lons_i = lons(ill);
                lats_i = lats(ill);
                % lower left origin
                
                [x_i,y_i] = geog2utm_nodisp(lons_i,lats_i,min(lons_i),min(lats_i)); 
            end
        end
%        [x,y] = geog2utm_nodisp(lons,lats,lon_min,lat_min); 
        % meshgrid after this
        
%         [x,y] = geog2utm_nodisp(lons,lats,lon_min,lat_min); % bottom left origin
%         x = reshape(x,size(lons));
%         y = reshape(y,size(lons));
        
        % UTM zones are 6 deg, do streamlines in 3x3 degs
        % utm conversion
        x = lons;
        y = lats;
        
        if res_x_mult == 1 % same res as vectors
            x_start = x;
            y_start = y;
            utest = u;
        else % higher res option
%             [x_start,y_start] = geog2utm_nodisp(lons_high,lats_high,lon_min,lat_min); % bottom left origin            
%             x_start = reshape(x_start,size(lons_high));
%             y_start = reshape(y_start,size(lons_high));
            x_start = lons_high;
            y_start = lats_high;
            utest = interp2(lons',lats',u',lons_high',lats_high','nearest')';
        end
        
        [row,col] = find(~isnan(utest)); % only run points that have data
        
        for k = 1:length(row)
            ik = row(k); % reversed indexing here?
            jk = col(k);
            startx = x_start(ik,jk);
            starty = y_start(ik,jk);
            try
                XY = stream2(x',y',u',v',startx,starty);
                % [stepsize, max_number_vertices]
                % If you do not specify a value, MATLAB® software uses the default:
                % Step size = 0.1 (one tenth of a cell)
                % Maximum number of vertices = 10000
%                 stepsize = 0.1;
%                 max_number_vertices = 10000;

                % limits to save computation time, doesn't change much
%                 stepsize = 0.1;
%                 max_number_vertices = 1000;
%                 XY = stream2(x',y',u',v',startx,starty,[stepsize, max_number_vertices]);
            catch % no streamline at that location
                continue
            end
            stline = XY{1};
            if length(stline(:)) < 20 % check to make sure it at least 10 points
                continue
            end
            if isnan(stline(10,1)) % only more than 10 points
                continue
            end
            
            % keep streamline data that is not nan
            stline_x = stline(:,1);
            stline_y = stline(:,2);
            last = find(isnan(stline_x),1);
            if ~isempty(last) % if no nans keep all the data
                stline_x = stline_x(1:last-1);
                stline_y = stline_y(1:last-1);
            end
            
            %%
            % winding angle calc start
            if plt_debug_streams == 1 && plt_debug == 1
                m_plot(stline_x,stline_y,'k')     % debug plot streamline
                drawnow                         % debug plot streamline
            end
            
            % initial angle
            eddy = 0; % not an eddy initially
            winding = 0; % winding angle
            dir = 0; % direction of winding
            baddir = 0; % initially no bad directions
            dx = stline_x(2) - stline_x(1);
            dy = stline_y(2) - stline_y(1);
            ang0 = atan2d(dy,dx);
            
            for ii = 2:length(stline_x)-1 % no angle at last data point
                dx = stline_x(ii+1) - stline_x(ii);
                dy = stline_y(ii+1) - stline_y(ii);
                ang1 = atan2d(dy,dx); % new angle
                ang_diff = ang1 - ang0; % winding angle is extra turning
                if dir == 0 % set with first point
                    dir = sign(ang_diff); % if returns zero will set next time
                end
                
                if ii > 2 % only after first winding angle
                    newdir = sign(ang_diff);
                    if newdir ~= dir % if winding changes direction
                        baddir = baddir + ang_diff;
                        if abs(baddir) > baddir_thres % break after this many degrees bad dirs
                            break
                        end
                    end
                    if newdir == dir % reset if goes back to correct direction before baddir_thres
                        baddir = 0;
                    end
                end
                
                if ang_diff > 180
                    ang_diff = ang_diff - 360;
                end
                if ang_diff < -180
                    ang_diff = ang_diff + 360;
                end
                
                ang0 = ang1; % reset angle
                winding = winding + ang_diff;
                
                if abs(winding) > winding_thres % if winding angle is close to 360
                    % check distance from starting point
                    dx = stline_x(ii) - stline_x(1);
                    dy = stline_y(ii) - stline_y(1);
                    d = sqrt(dx^2+dy^2);
                    if eddy == 1 % if already eddy see if distance is closer
                        if d > d1 % if distance increases break and save
                            break
                        end
                    end % otherwise distance decreased and keep going
                    if d < d_thres % now it is an eddy!
                        eddy = 1;
                        d1 = d; % save this distance
                    end % we should still keep going to see if d decreases
                end % otherwise keep going
            end
            
            if eddy == 1
                i_streams = i_streams + 1;
                streamsx{i_streams} = stline_x(1:ii); % up to last pt before break
                streamsy{i_streams} = stline_y(1:ii);
                streamsx_center{i_streams} = mean(streamsx{i_streams});
                streamsy_center{i_streams} = mean(streamsy{i_streams});
                winding_dir(i_streams) = dir; % math notation
                
                % angular velocity in deg/sec
                u1 = interp2(x',y',u',streamsx{i_streams},streamsy{i_streams});
                v1 = interp2(x',y',v',streamsx{i_streams},streamsy{i_streams});
                vv = sqrt(u1.^2 + v1.^2); % in m/s
                vv = vv(1:end-1); % last velocity doesn't add to streamline
                dx = diff(streamsx{i_streams});
                dy = diff(streamsy{i_streams});
                dr = sqrt(dx.^2 + dy.^2); % in km
                dt = (dr./vv)*1000; % time in seconds
                dt = sum(dt);
                angular_vel(i_streams) = winding/dt; % degs/second
%               
                % debug plot
                if plt_debug == 1
                    m_plot(streamsx{i_streams},streamsy{i_streams},...
                        'm','linewidth',2)
                    m_plot(streamsx_center{i_streams},streamsy_center{i_streams},...
                        'rs','linewidth',2)
%                     drawnow
                end
%                 disp('found eddy') % debug console info
                
            end
        end % each streamline
        
        %%
        if plt_debug == 1
            for iplt = 1:i_streams
                m_plot(streamsx{iplt},streamsy{iplt},...
                    'm','linewidth',2)
                m_plot(streamsx_center{iplt},streamsy_center{iplt},...
                    'rs','linewidth',2)
                %         drawnow
            end
        end
        %% streamline clustering
        eddy_nums = 1:i_streams; % eddy numbers
        for i1 = 1:i_streams - 1
            if eddy_nums(i1) ~= i1
                continue
            end
            if plt_debug == 1
                
                m_plot(streamsx{i1},streamsy{i1},...
                    'm','linewidth',2)
                m_plot(streamsx_center{i1},streamsy_center{i1},...
                    'rs','linewidth',2)
                %         drawnow
            end
            
            % check to see if centers are close enough
            for j1 = i1:i_streams
                if i1 == j1
                    continue
                end
                
                dx = streamsx_center{i1} - streamsx_center{j1};
                dy = streamsy_center{i1} - streamsy_center{j1};
                d = sqrt(dx^2 + dy^2);
                
                % if centers are close and winding direction is the same
                if d < param_center_dist && (winding_dir(i1) == winding_dir(j1)) % then eddy are in the same cluster
                    eddy_nums(j1) = i1; % same eddy as i
%                     disp('same eddy') % debug console info
                end
                
            end
        end
        
        %% s
        i1 = unique(eddy_nums);
        eddies = 0;
        for j1 = 1:length(i1)
            k = find(eddy_nums == i1(j1));
            eddies = eddies + 1;
            ptsx = [];
            ptsy = [];
            eddy_angular_vel1 = [];
            for ii = 1:length(k)
                ptsx = [ptsx ; streamsx{k(ii)}];
                ptsy = [ptsy ; streamsy{k(ii)}];
                eddy_angular_vel1 = [eddy_angular_vel1 angular_vel(k(ii))];
            end
            eddy_center_x(eddies) = mean(ptsx);
            eddy_center_y(eddies) = mean(ptsy);
            eddy_dir(eddies) = winding_dir(k(1)); % winding dir, they are all the same
            eddy_angular_vel(eddies) = mean(eddy_angular_vel1); % avg angular velocity 
            eddy_streamlines(eddies) = length(k);
            %%
            %     [d,l] = ellipse_fit(ptsx,ptsy); % d is directions, l is lengths
            % these three lines are the function
            M = cov(ptsx,ptsy);
            [d,l] = eig(M); % columns are eig vecs
            l = diag(l);
            
            dx1 = d(1,1);
            dy1 = d(2,1);
            l1 = sqrt(l(1)); % sqrt?
            l2 = sqrt(l(2)); % sqrt?
             
            ellipse_theta = atan2d(dy1,dx1);
            eddy_ellipse_theta(eddies) = ellipse_theta;
            eddy_length_x(eddies) = l1;
            eddy_length_y(eddies) = l2;
            %     d(:,1);
            %     d(:,2);
            if plt_debug == 1 % ellipse plots
                it = 1:360;
                ellipse_x = l1*cosd(it);
                ellipse_y = l2*sind(it);
                % Create rotation matrix
                R = [cosd(ellipse_theta) -sind(ellipse_theta); sind(ellipse_theta) cosd(ellipse_theta)];
                % Rotate points
                eddypts = R*[ellipse_x ; ellipse_y];
                ellipse_x = eddypts(1,:);
                ellipse_y = eddypts(2,:);
            
                m_plot(eddy_center_x(eddies),eddy_center_y(eddies),...
                    'kx','linewidth',4,'markersize',12)
                m_plot(eddy_center_x(eddies) + ellipse_x,...
                    eddy_center_y(eddies) + ellipse_y,'b.')
            end
        end
        
        % show time of calculation
        numfiles_i = numfiles_i + 1; % file count
        toc
        
        if eddies == 0
            eddy_center_x = [];
            eddy_center_y = [];
            eddy_dir = [];
            eddy_angular_vel = [];
            eddy_streamlines = [];
            eddy_length_x = [];
            eddy_length_y = [];
            eddy_ellipse_theta = [];
        end
        % save data
        save(snmat,'eddy_center_x','eddy_center_y',...
            'eddy_dir','eddy_angular_vel','eddy_streamlines',...
            'eddy_length_x','eddy_length_y','eddy_ellipse_theta');
        
        if plt_debug == 1 % save plot
            savefig(sn)
            print(sn,'-dpng')
            if numfiles ~= 1
                close(h)
            end
        end
        
        disp([num2str(eddies) ' eddies found'])
        for i_eddy = 1:eddies
            disp(['eddy #' num2str(i_eddy) ' has ' ...
                num2str(eddy_streamlines(i_eddy)) ' streamlines'])
        end
        disp(' ')
        
        
        
    end
    if numfiles > 0
        if numfiles_i > numfiles
            return
        end
    end
end




