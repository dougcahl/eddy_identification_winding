%% function [] = eddy_subroutine(data,params)
function [] = eddy_subroutine(params,u,v,time,lons,lats,lon,lat)
%
% Douglas Cahl, PhD 
% University of South Carolina 2023


% parameters
save_data_name = params.save_data_name;
save_fig_name = params.save_fig_name;
save_streams = params.save_streams;
dij         = params.dij; % streamlines in areas of lons(i-dij:i+dij,j-dij:j+dij)
% for 1km resolution, dij = 100 for max diameter streamline of ~200 km

lonlatlimits = params.lonlatlimits; % (=1) enable limiting to subregion, set 0 to do entire region
if lonlatlimits
    lon_min = params.lon_min; 
    lon_max = params.lon_max;
    lat_min = params.lat_min;
    lat_max = params.lat_max;
end

plt_debug = params.plt_debug;                   % debug plots
plt_debug_streams = params.plt_debug_streams;   % streamline plots
plt_quiver = params.plt_quiver;                 % plot quiver
plt_debug_sv = params.plt_debug_sv;             % save debug plot
highrescoast = params.highrescoast;             % high resolution coastline
bathy = params.bathy;                           % plot bathy levels
levels = params.levels;                         % isobaths 
plt_close = params.plt_close;                   % close figures after running
print_res = params.print_res;                   % resolution of plot 
plt_each_eddy = params.plt_each_eddy;           % zooms and plots figure for each eddy
plt_min_streams = params.plt_min_streams;
dlon = params.radar_km_resolution/6;
dlat = dlon;

winding_thres = params.winding_thres;           % how much winding (in degrees) to be a eddy
d_thres = params.d_thres;                       % streamline starting and ending max distance in km
baddir_thres = params.baddir_thres;             % breaks at x degrees of bad dirs
param_center_dist = params.param_center_dist;   % how many km are clusters together
new_dist_thres = params.new_dist_thres;         % Cahl et al., 2023 method of cluster distance


if lonlatlimits
    
    j = lons >= lon_min & lons <= lon_max & lats >= lat_min & lats <= lat_max;
    u = u(j);
    v = v(j);
    lons = lons(j);
    lats = lats(j);
    
    j = lon >= lon_min & lon <= lon_max;
    lon1 = lon(j);
    j = lat >= lat_min & lat <= lat_max;
    lat1 = lat(j);
    [lats1,~] = meshgrid(lat1,lon1);
    lons = reshape(lons,size(lats1));
    lats = reshape(lats,size(lats1));
    u = reshape(u,size(lats1));
    v = reshape(v,size(lats1));
    
else
    lon_min = double(floor(10*min(lons(:)))/10);
    lon_max = double(ceil(10*max(lons(:)))/10);
    lat_min = double(floor(10*min(lats(:)))/10);
    lat_max = double(ceil(10*max(lats(:)))/10);
end

%%
if plt_debug_streams
    h = figure('units','normalized','outerposition',[0.05 0.05 .9 .9]);
    m_proj('Mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max])
    m_grid();
    hold on
    
    if plt_quiver
        %     m_quiver(lons,lats,u,v)
        U = sqrt(u.^2 + v.^2);
        m_pcolor(lons,lats,U)
        m_quiver(lons,lats,u./U,v./U,0.5,'k')
        c = colorbar;
        ylabel(c,'Speed (m/s)')
        c1 = clim;
    end
  
    if bathy
        m_etopo2('contour',levels,'edgecolor','k','linewidth',0.5,'linestyle','-','ShowText','on');
        if plt_quiver
            clim(c1);
        end
    end
    
    if highrescoast
        m_gshhs_h('patch',[0.8 0.8 0.8]);
    else
        m_gshhs_l('patch',[0.8 0.8 0.8]);
    end
end

%%

i_streams = 0;
% UTM zones are 6 deg, should do streamlines in smaller <3x3 degs

[row,col] = find(~isnan(u)); % only run points that have data

% dij = 100;
for k = 1:length(row)
    ik = row(k); % reversed indexing here?
    jk = col(k);

    start_lon = lons(ik,jk);
    start_lat = lats(ik,jk);

    ik = max([1,ik-dij]) : min([size(lons,1),ik+dij]);
    jk = max([1,jk-dij]) : min([size(lons,2),jk+dij]);

    lons1 = lons(ik,jk);
    lats1 = lats(ik,jk);
    
    u1 = double(u(ik,jk));
    v1 = double(v(ik,jk));
    
    % UTM with starting position at 0,0
%     [x1,y1,UTMzone]=geog2utm_nodisp(lons1,lats1,start_lon,start_lat); % in km
    [x1,y1,UTMzone]=geog2utm_nodisp(lons1,lats1); % in km
    [startx,starty,UTMzone]=geog2utm_nodisp(start_lon,start_lat,UTMzone); % in km
    x1 = reshape(x1,size(lons1));
    y1 = reshape(y1,size(lons1));

%     tic
%     dx = diff(x1);
%     dx = mean(dx(:));
%     dy = diff(y1');
%     dy = mean(dy(:));
%     ddx = min(dx,dy);
%     [y2,x2] = meshgrid(min(y1(:)):ddx:max(y1(:)),min(x1(:)):ddx:max(x1(:)));
%     u1 = griddata(x1,y1,u1,x2,y2);
%     v1 = griddata(x1,y1,v1,x2,y2);
% %     f = scatteredInterpolant(x1(:),y1(:),double(u1(:)),'linear','none');
% %     u1 = f(x2,y2);
% %     f = scatteredInterpolant(x1(:),y1(:),double(v1(:)));
% %     v1 = f(x2,y2);
%     x1 = x2;
%     y1 = y2;
%     clear x2 y2 f
%     toc
    



%     % starting position at 0,0 
%     startx = 0;
%     starty = 0;
    
    try
%         XY = stream2(x1',y1',u1',v1',startx,starty);
%         [X,Y] = stream2_dc(x1,y1,u1,v1,startx,starty);
        [X,Y] = stream2_dc(x1,y1,u1,v1,startx,starty,.2,2000);
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
%     stline = XY{1};
%     if length(stline(:)) < 20 % check to make sure it at least 10 points
%         continue
%     end
%     if isnan(stline(10,1)) % only more than 10 points
%         continue
%     end

    if length(X) < 20 % check to make sure it at least 10 points
        continue
    end
    if isnan(X(10)) % only more than 10 points
        continue
    end
    
    % keep streamline data that is not nan
%     stline_x = stline(:,1);
%     stline_y = stline(:,2);
%     last = find(isnan(stline_x),1);
%     if ~isempty(last) % if no nans keep all the data
%         stline_x = stline_x(1:last-1);
%         stline_y = stline_y(1:last-1);
%     end

    last = find(isnan(X),1);
    if ~isempty(last) % if no nans keep all the data
        stline_x = X(1:last-1);
        stline_y = Y(1:last-1);
    else
        stline_x = X;
        stline_y = Y;
    end

    
    %%
    % winding angle calc start
    if plt_debug_streams && plt_debug
        [stline_lat,stline_lon]=utm2ll(stline_x*1000,stline_y*1000,UTMzone);
        m_plot(stline_lon,stline_lat,'k')   % debug plot streamline
        drawnow                         % debug plot streamline
    end
    
    % initial angle
    eddy = 0;       % not an eddy initially
    winding = 0;    % winding angle
    dir = 0;        % direction of winding
    baddir = 0;     % initially no bad directions
    dx = stline_x(2) - stline_x(1);
    dy = stline_y(2) - stline_y(1);
    ang0 = atan2d(dy,dx);
    
    for ii = 2:length(stline_x)-1 % no angle at last data point
        dx = stline_x(ii+1) - stline_x(ii);
        dy = stline_y(ii+1) - stline_y(ii);
        ang1 = atan2d(dy,dx);       % new angle
        ang_diff = ang1 - ang0;     % winding angle is extra turning
        if ang_diff > 180
            ang_diff = ang_diff - 360;
        end
        if ang_diff < -180
            ang_diff = ang_diff + 360;
        end
        if dir == 0                 % set with first point
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
        
        ang0 = ang1; % reset angle
        winding = winding + ang_diff;
        
        if abs(winding) > winding_thres % if winding angle is close to 360
            % check distance from starting point
            dx = stline_x(ii) - stline_x(1);
            dy = stline_y(ii) - stline_y(1);
            d = sqrt(dx^2+dy^2);
            if eddy % if already eddy see if distance is closer
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
    
    if eddy
        i_streams = i_streams + 1;
        streamsx{i_streams} = stline_x(1:ii); % up to last pt before break
        streamsy{i_streams} = stline_y(1:ii);
        streamsx_center{i_streams} = mean(streamsx{i_streams});
        streamsy_center{i_streams} = mean(streamsy{i_streams});
        [streamslat{i_streams},streamslon{i_streams}]= utm2ll(streamsx{i_streams}*1000 , streamsy{i_streams}*1000,UTMzone);
        [streamslat_center{i_streams},streamslon_center{i_streams}]= utm2ll(streamsx_center{i_streams}*1000 , streamsy_center{i_streams}*1000,UTMzone);



        winding_dir(i_streams) = dir; % math notation

        % angular velocity in deg/sec
        f11u = scatteredInterpolant(x1(:),y1(:),u1(:));
        f11v = scatteredInterpolant(x1(:),y1(:),v1(:));
        u2 = f11u(streamsx{i_streams},streamsy{i_streams});
        v2 = f11v(streamsx{i_streams},streamsy{i_streams});
        vv = sqrt(u2.^2 + v2.^2);   % in m/s
        vv = vv(1:end-1);           % last velocity doesn't add to streamline
        dx = diff(streamsx{i_streams});
        dy = diff(streamsy{i_streams});
        dr = sqrt(dx.^2 + dy.^2);   % in km
        dt = (dr./vv)*1000;         % time in seconds
        dt = sum(dt);
        angular_vel(i_streams) = winding/dt; % degs/second

        

        
    end
end % each streamline

%% streamline clustering
eddy_nums = 1:i_streams; % eddy numbers
for i1 = 1:i_streams - 1
    if eddy_nums(i1) ~= i1
        continue
    end

    % check to see if centers are close enough
    for j1 = i1:i_streams
        if i1 == j1
            continue
        end

        if new_dist_thres % Cahl et al., 2023 method (if center is within another streamline, they are the same eddy)
            in1 = inpolygon(streamslon_center{j1},streamslat_center{j1},...
                streamslon{i1},streamslat{i1});
            in2 = inpolygon(streamslon_center{i1},streamslat_center{i1},...
                streamslon{j1},streamslat{j1});
            in11 = inpolygon(streamslon{j1},streamslat{j1},...
                streamslon{i1},streamslat{i1});
            in22 = inpolygon(streamslon{i1},streamslat{i1},...
                streamslon{j1},streamslat{j1});
            if (in1 || in2 || sum(in11) > 1 || sum(in22) > 1) % && (winding_dir(i1) == winding_dir(j1)) % if center is within another streamline, they are the same eddy
                eddy_nums(j1) = i1; % same eddy as i
            else
                [dx,dy] = geog2utm_nodisp(streamslon_center{i1},streamslat_center{i1},streamslon_center{j1},streamslat_center{j1});
                d = sqrt(dx^2 + dy^2); % in km
                % if centers are close and winding direction is the same
                if d < param_center_dist % && (winding_dir(i1) == winding_dir(j1)) % then same eddy/cluster
                    eddy_nums(j1) = i1; % same eddy as i
                end
            end
        
        else % Sadarjoen and Post 2000 method
            [dx,dy] = geog2utm_nodisp(streamslon_center{i1},streamslat_center{i1},streamslon_center{j1},streamslat_center{j1});
            d = sqrt(dx^2 + dy^2); % in km
            % if centers are close and winding direction is the same
            if d < param_center_dist && (winding_dir(i1) == winding_dir(j1)) % then same eddy/cluster
                eddy_nums(j1) = i1; % same eddy as i
            end
        end

    end % for each other streamline
end % for each streamline


%%
if plt_debug && ~plt_each_eddy
    h = figure('units','normalized','outerposition',[0.05 0.05 .9 .9]);
    m_proj('Mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max])
    m_grid();
    hold on
    if plt_quiver
        %     m_quiver(lons,lats,u,v)
        U = sqrt(u.^2 + v.^2);
        m_pcolor(lons,lats,U)
        m_quiver(lons,lats,u./U,v./U,0.5,'k')
        c = colorbar;
        ylabel(c,'Speed (m/s)')
        c1 = clim;
    end
    if bathy
        m_etopo2('contour',levels,'edgecolor','k','linewidth',0.5,'linestyle','-','ShowText','on');
        if plt_quiver
            clim(c1);
        end
    end
    if highrescoast
        m_gshhs_h('patch',[0.8 0.8 0.8]);
    else
        m_gshhs_l('patch',[0.8 0.8 0.8]);
    end
end
%% calculate eddy parameters
i1 = unique(eddy_nums);
eddies = 0;

eddy_streamlines_lat = [];
eddy_streamlines_lon = [];
for j1 = 1:length(i1)
    k = find(eddy_nums == i1(j1));
    eddies = eddies + 1;
    ptsx = [];
    ptsy = [];
    eddy_angular_vel1 = [];
    for ii = 1:length(k)
        if save_streams
            eddy_streamlines_lon{eddies,ii} = streamslon{k(ii)};
            eddy_streamlines_lat{eddies,ii} = streamslat{k(ii)};
        end

        ptsx = [ptsx ; streamslon{k(ii)}];
        ptsy = [ptsy ; streamslat{k(ii)}];
        eddy_angular_vel1 = [eddy_angular_vel1 angular_vel(k(ii))];
    end
    [ptsx,ptsy,UTMzone] = geog2utm_nodisp(ptsx,ptsy);
    eddy_center_x(eddies) = mean(ptsx);
    eddy_center_y(eddies) = mean(ptsy);
    [eddy_center_lat1,eddy_center_lon1]=utm2ll(eddy_center_x(eddies)*1000,...
            eddy_center_y(eddies)*1000,UTMzone);
    eddy_center_lon(eddies) = eddy_center_lon1;
    eddy_center_lat(eddies) = eddy_center_lat1;

    eddy_dir(eddies) = winding_dir(k(1)); % winding dir, they are all the same
    eddy_angular_vel(eddies) = mean(eddy_angular_vel1); % avg angular velocity
    eddy_streamlines(eddies) = length(k);
    %
    M = cov(ptsx,ptsy);
    [d,l] = eig(M); % columns are eig vecs
    l = diag(l);

    dx1 = d(1,1);
    dy1 = d(2,1);
    l1 = sqrt(2*l(1));
    l2 = sqrt(2*l(2));

    ellipse_theta = atan2d(dy1,dx1);
    eddy_ellipse_theta(eddies) = ellipse_theta;
    eddy_length_x(eddies) = l1;
    eddy_length_y(eddies) = l2;
  
    % streamline and ellipse plots
    if plt_debug && eddy_streamlines(eddies) > plt_min_streams 
        if plt_each_eddy % new figure
            h = figure('units','normalized','outerposition',[0.05 0.05 .9 .9]);
            m_proj('Mercator','lon',[eddy_center_lon1-dlon eddy_center_lon1+dlon],...
                'lat',[eddy_center_lat1-dlat eddy_center_lat1+dlat])
            m_grid();
            hold on
            if plt_quiver
                %     m_quiver(lons,lats,u,v)
                U = sqrt(u.^2 + v.^2);
                m_pcolor(lons,lats,U)
                m_quiver(lons,lats,u./U,v./U,0.5,'k')
                c = colorbar;
                ylabel(c,'Speed (m/s)')
                c1 = clim;
            end
            if bathy
                m_etopo2('contour',levels,'edgecolor','k','linewidth',0.5,'linestyle','-','ShowText','on');
                if plt_quiver
                    clim(c1);
                end
            end
            if highrescoast
                m_gshhs_h('patch',[0.8 0.8 0.8]);
            else
                m_gshhs_l('patch',[0.8 0.8 0.8]);
            end
        end

        for ii = 1:length(k)
            m_plot(streamslon{k(ii)},streamslat{k(ii)},'m','linewidth',0.25)
            m_plot(streamslon_center{k(ii)},streamslat_center{k(ii)},'rx','linewidth',1)
        end

        it = 1:360;
        ellipse_x = l1*cosd(it);
        ellipse_y = l2*sind(it);
        % Create rotation matrix
        R = [cosd(ellipse_theta) -sind(ellipse_theta); sind(ellipse_theta) cosd(ellipse_theta)];
        % Rotate points
        eddypts = R*[ellipse_x ; ellipse_y];
        ellipse_x = eddypts(1,:);
        ellipse_y = eddypts(2,:);
        
        m_plot(eddy_center_lon1,eddy_center_lat1,...
            'kx','linewidth',2,'markersize',12)
        [eddy_center_x,eddy_center_y,UTMzone] = geog2utm_nodisp(eddy_center_lon1,eddy_center_lat1);
        [eddy_center_lats_ellipse,eddy_center_lons_ellipse]=utm2ll((eddy_center_x + ellipse_x)*1000,...
            (eddy_center_y + ellipse_y)*1000,UTMzone);
        m_plot([eddy_center_lons_ellipse ; eddy_center_lons_ellipse(1)],[eddy_center_lats_ellipse ; eddy_center_lats_ellipse(1)],...
            'color',[0.8 0.8 0.8],'linewidth',1)
        drawnow

        if plt_debug_sv && plt_debug && plt_each_eddy
        %     savefig(save_fig_name)
            print([save_fig_name '_eddy' num2str(j1)],'-dpng',['-r' num2str(print_res)])
            if plt_close
                close(h)
            end
        end

    end
end

%% save data
if eddies == 0
    eddy_center_lon = [];
    eddy_center_lat = [];
    eddy_dir = [];
    eddy_angular_vel = [];
    eddy_streamlines = [];
    eddy_length_x = [];
    eddy_length_y = [];
    eddy_ellipse_theta = [];
end

save(save_data_name,'eddy_center_lat','eddy_center_lon',...
    'eddy_dir','eddy_angular_vel','eddy_streamlines',...
    'eddy_length_x','eddy_length_y','eddy_ellipse_theta',...
    'eddy_streamlines_lat','eddy_streamlines_lon',...
    'UTMzone','time');

%% save plot
if plt_debug_sv && plt_debug && ~plt_each_eddy
%     savefig(save_fig_name)
    print(save_fig_name,'-dpng',['-r' num2str(print_res)])
    if plt_close
        close(h)
    end
end

disp([num2str(eddies) ' eddies found'])
for i_eddy = 1:eddies
    disp(['eddy #' num2str(i_eddy) ' has ' ...
        num2str(eddy_streamlines(i_eddy)) ' streamlines'])
end
disp(' ')
