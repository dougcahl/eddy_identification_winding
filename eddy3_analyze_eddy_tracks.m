%% analyze eddy tracks
%
% Douglas Cahl, PhD 
% University of South Carolina 2023

%% parameters
params_eddy3;

%% directory setup
name_dir = tracks_dir;     % directory of data
name_pre  = [region '_' num2str(radar_km_resolution) 'km_']; % name prefix, e.g., data2_1.mat , data2_2.mat , ...
tracksfile = [name_dir name_pre 'tracks'];
if ~exist([tracksfile '.mat'],'file')
    return
end
load(tracksfile)



%%
disp(['max gap between eddy identification in timesteps = ' num2str(eddy_track_time_param)])% - 1)])
disp(['number of tracks = ' num2str(length(timegap))])
for i = 1:length(timegap)
    eddy_track_length(i) = length(direction{i});
    intensity(i) = mean(num_streams{i});
    t(i) = mean(Time{i});
    d = eig1{i} + eig2{i};
    diameter(i) = mean(d); % 2x average of major and minor axes
    maxdiameter(i) = max(d); % 2x max of major and minor axes
    mindiameter(i) = min(d); % 2x max of major and minor axes

    if eddy_track_length(i) > 1
        lon = lon_center{i};
        lat = lat_center{i};
        [x,y,UTMzone] = geog2utm_nodisp(lon,lat);
        
        dx = x(end)-x(1);
        dy = y(end)-y(1);
        dist1 = sqrt(dx^2+dy^2)*1000;
        dt = Time{i}(end)-Time{i}(1);
        dt = dt*24*3600;
        dist(i) = dist1;
        speed(i) = dist1/dt;
    else
        dist(i) = nan;
        speed(i) = nan;
    end

end

disp(['max eddy track length in timesteps = ' num2str(max(eddy_track_length))])


%% tracks by length plots

edl = eddy_track_length;


[edl,b] = sort(edl,'descend');
counter = 0;
w3 = [];
for i = 1:length(edl)
    j = b(i);
    
    if eddy_track_length(j) < min_steps 
        continue
    end
    counter = counter + 1;
    
    w1 = mean(omega{j},'omitnan');
    w2 = w1*3600*24/360;
    w3(i) = w2;
    
    lon = lon_center{j};
    lat = lat_center{j};
    
    if plt
        h = figure('color','w');
        m_proj('Mercator','lon',[mean(lon)-dlon mean(lon)+dlon],...
                    'lat',[mean(lat)-dlat mean(lat)+dlat])
      
        if highrescoast
            m_gshhs_h('patch',[0.8 0.8 0.8]); % coastline
        else
            m_gshhs_l('patch',[0.8 0.8 0.8]); % coastline
        end
        levels = [-10 -25 -50 -100 -250 -500 -1000 -2500];
        m_etopo2('contour',levels,'edgecolor',0.5*[1 1 1],...
            'linewidth',0.5,'linestyle','-','ShowText','off');
        m_grid(); 
        hold on
        
        if w1 > 0
            m_plot(lon,lat,'-b.') % CCW (Northern: cyclonic/upwelling ; Southern, anticyclonic/downwelling)
        else
            m_plot(lon,lat,'-r.') % CW (Northern: anticyclonic/downwelling ; Southern, cyclonic/upwelling)
        end
        
        m_plot(lon(1),lat(1),'o','LineWidth',2,'color',[0 .8 0])
        m_plot(lon(end),lat(end),'kx','LineWidth',2)
        
        title(['track lenth = ' num2str(eddy_track_length(j)) ' timesteps  ;   ' ...
            num2str(round(w2,1)) ' rot/day'])
        
        print([tracksfile '_' num2str(counter)],'-dpng',['-r' num2str(print_res)])
        if plt_close
            close(h)
        end
    end
end



%% histogram plots
ed1 = eddy_track_length*radar_dt; % conver to hours

if plt
    h = figure('color','w');
    subplot(231)
    histogram(ed1)
    xlabel('Eddy duration (hours)')
    ylabel('Number of eddies')
    
    
    subplot(232)
    histogram(intensity)
    xlabel('intensity (avg streamlines)')
    ylabel('Number of eddies')
    box off
    grid on
    
    subplot(233)
    histogram(speed)
    xlabel('speed (m/s)')
    box off
    grid on
    
    subplot(234)
    histogram(diameter)
    xlabel('avg diameter (km)')
    box off
    grid on
    
    subplot(235)
    histogram(mindiameter)
    xlabel('min diameter (km)')
    box off
    grid on
    
    subplot(236)
    histogram(maxdiameter)
    xlabel('max diameter (km)')
    box off
    grid on
    
    print([tracksfile '_hist'],'-dpng',['-r' num2str(print_res)])
    if plt_close
        close(h)
    end
end




