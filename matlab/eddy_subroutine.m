%% function [] = eddy_subroutine(data,params)
function [] = eddy_subroutine(params,u,v,time,lons,lats,lon,lat)
%
% Douglas Cahl, PhD 
% University of South Carolina 2023


% parameters
save_name = params.save_name;

lonlatlimits = params.lonlatlimits; % (=1) enable limiting to subregion, set 0 to do entire region
if lonlatlimits
    lon_min = params.lon_min; 
    lon_max = params.lon_max;
    lat_min = params.lat_min;
    lat_max = params.lat_max;
end

plt_debug = params.plt_debug;                   % debug plots
plt_debug_streams = params.plt_debug_streams;   % streamline plots
plt_debug_sv = params.plt_debug_sv;             % save debug plot
highrescoast = params.highrescoast;             % high resolution coastline
bathy = params.bathy;                           % plot bathy levels
levels = params.levels;                         % isobaths 
plt_close = params.plt_close;                   % close figures after running
print_res = params.print_res;                   % resolution of plot 

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
if plt_debug
    h = figure('units','normalized','outerposition',[0.05 0.05 .9 .9]);
    m_proj('Mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max])
    m_grid();
    hold on

    %     m_quiver(lons,lats,u,v)
    U = sqrt(u.^2 + v.^2);
    m_pcolor(lons,lats,U)
    m_quiver(lons,lats,u./U,v./U,0.5,'k')
    c = colorbar;
    ylabel(c,'Speed (m/s)')
    c1 = clim;

    if bathy
        m_etopo2('contour',levels,'edgecolor','k','linewidth',0.5,'linestyle','-','ShowText','on');
        clim(c1);
    end
    
    if highrescoast
        m_gshhs_h('patch',[0.8 0.8 0.8]);
    else
        m_gshhs_l('patch',[0.8 0.8 0.8]);
    end
end

%%

i_streams = 0;
% UTM zones are 6 deg, should do streamlines in 3x3 degs
[x,y,UTMzone]=geog2utm_nodisp(lons,lats); % in km
x = reshape(x,size(lons));
y = reshape(y,size(lons));

x_start = x;
y_start = y;
utest = u;

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
    if plt_debug_streams && plt_debug
        m_plot(stline_x,stline_y,'k')   % debug plot streamline
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
        winding_dir(i_streams) = dir; % math notation

        % angular velocity in deg/sec
        f11u = scatteredInterpolant(x(:),y(:),u(:));
        f11v = scatteredInterpolant(x(:),y(:),v(:));
        u1 = f11u(streamsx{i_streams},streamsy{i_streams});
        v1 = f11v(streamsx{i_streams},streamsy{i_streams});
        vv = sqrt(u1.^2 + v1.^2);   % in m/s
        vv = vv(1:end-1);           % last velocity doesn't add to streamline
        dx = diff(streamsx{i_streams});
        dy = diff(streamsy{i_streams});
        dr = sqrt(dx.^2 + dy.^2);   % in km
        dt = (dr./vv)*1000;         % time in seconds
        dt = sum(dt);
        angular_vel(i_streams) = winding/dt; % degs/second
        
    end
end % each streamline

%% plot eddy streamlines 
if plt_debug
    for iplt = 1:i_streams
        [streams_lat,streams_lon]=utm2ll(streamsx{iplt}*1000,streamsy{iplt}*1000,UTMzone);
        m_plot(streams_lon,streams_lat,'m','linewidth',0.25)
        
        [streamsy_center_lat,streamsx_center_lon]=utm2ll(streamsx_center{iplt}*1000,streamsy_center{iplt}*1000,UTMzone);
        m_plot(streamsx_center_lon,streamsy_center_lat,'rx','linewidth',1)
    end
end
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
            in1 = inpolygon(streamsx_center{j1},streamsy_center{j1},...
                streamsx{i1},streamsy{i1});
            in2 = inpolygon(streamsx_center{i1},streamsy_center{i1},...
                streamsx{j1},streamsy{j1});
            if in1 || in2 % if center is within another streamline, they are the same eddy
                eddy_nums(j1) = i1; % same eddy as i
            end

        else % Sadarjoen and Post 2000 method
            dx = streamsx_center{i1} - streamsx_center{j1};
            dy = streamsy_center{i1} - streamsy_center{j1};
            d = sqrt(dx^2 + dy^2); % in km
            % if centers are close and winding direction is the same
            if d < param_center_dist && (winding_dir(i1) == winding_dir(j1)) % then same eddy/cluster
                eddy_nums(j1) = i1; % same eddy as i
            end
        end

    end % for each other streamline
end % for each streamline

%% second round: Cahl et al., 2023 method (if center is within another streamline, they are the same eddy)
if new_dist_thres 
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
            
        end
        eddy_center_x(eddies) = mean(ptsx);
        eddy_center_y(eddies) = mean(ptsy);
        %
        M = cov(ptsx,ptsy);
        [d,l] = eig(M); % columns are eig vecs
        l = diag(l);
    
        dx1 = d(1,1);
        dy1 = d(2,1);
        l1 = sqrt(2*l(1));
        l2 = sqrt(2*l(2));
    
        ellipse_theta = atan2d(dy1,dx1);
      
        it = 0:15:360;
        ellipse_x = l1*cosd(it);
        ellipse_y = l2*sind(it);
        % Create rotation matrix
        R = [cosd(ellipse_theta) -sind(ellipse_theta); sind(ellipse_theta) cosd(ellipse_theta)];
        % Rotate points
        eddypts = R*[ellipse_x ; ellipse_y];
        ellipse_x = eddypts(1,:);
        ellipse_y = eddypts(2,:);
        poly_x =  mean(ptsx) + ellipse_x;
        poly_y =  mean(ptsy) + ellipse_y;
    
        for j2 = 1:length(i1)
            if i1(j1) == i1(j2)
                continue
            end
            k2 = find(eddy_nums == i1(j2));
            for k3 = 1:length(k2)
                in = inpolygon(streamsx_center{k2(k3)},streamsy_center{k2(k3)},...
                    poly_x,poly_y);
                if in % if center is within another streamline, they are the same eddy
                    eddy_nums(k2(k3)) = i1(j1); % same eddy as i
                end
            end
        end
        
    end
end

%% cluster and calculate eddy parameters
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
  
    if plt_debug % ellipse plots
        it = 1:360;
        ellipse_x = l1*cosd(it);
        ellipse_y = l2*sind(it);
        % Create rotation matrix
        R = [cosd(ellipse_theta) -sind(ellipse_theta); sind(ellipse_theta) cosd(ellipse_theta)];
        % Rotate points
        eddypts = R*[ellipse_x ; ellipse_y];
        ellipse_x = eddypts(1,:);
        ellipse_y = eddypts(2,:);
        
        [eddy_center_lat,eddy_center_lon]=utm2ll(eddy_center_x(eddies)*1000,...
            eddy_center_y(eddies)*1000,UTMzone);
        m_plot(eddy_center_lon,eddy_center_lat,...
            'kx','linewidth',2,'markersize',12)
        [eddy_center_lats,eddy_center_lons]=utm2ll((eddy_center_x(eddies) + ellipse_x)*1000,...
            (eddy_center_y(eddies) + ellipse_y)*1000,UTMzone);
        m_plot([eddy_center_lons ; eddy_center_lons(1)],[eddy_center_lats ; eddy_center_lats(1)],...
            'color',[0.8 0.8 0.8],'linewidth',1)
        drawnow
    end
end

%% save data
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

save(save_name,'eddy_center_x','eddy_center_y',...
    'eddy_dir','eddy_angular_vel','eddy_streamlines',...
    'eddy_length_x','eddy_length_y','eddy_ellipse_theta','UTMzone','time');

%% save plot
if plt_debug_sv && plt_debug
%     savefig(save_name)
    print(save_name,'-dpng',['-r' num2str(print_res)])
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
