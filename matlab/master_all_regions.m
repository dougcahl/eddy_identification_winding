%% master_all_regions
% calls master1(region_num,radar_km_resolution,avg_25hr)
% where 
% region_num = 1    region = 'USEGC'    % 1,2,6 km                              fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'
% region_num = 2    region = 'GAK'      % 1,2,6 km                              fullregionname = 'HFRADAR_US_Gulf_of_Alaska'
% region_num = 3    region = 'AKNS'     % 6 km                                  fullregionname = 'HFRADAR_Alaska_-_North_Slope'
% region_num = 4    region = 'USWC'     % 500m 1,2,6 km (500m is not supported)    fullregionname = 'HFRADAR_US_West_Coast'
% region_num = 5    region = 'USHI'     % 1,2,6 km                              fullregionname = 'HFRADAR_US_Hawaii'
% region_num = 6    region = 'PRVI'     % 2,6 km             fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands'
% example: pythong3 master1.py 1 6 1

% run with:   ./eddy_all_regions.sh > master.log
% find with:  sudo ps -af | grep master1
% kill with:  sudo pkill -f 'master1.py'
% check logs: ls *.txt | xargs -I{} tail -1 {}

%     matrun = '!matlab -r "';
matrun = '!matlab -nosplash -nodesktop -r "';
matvars = ['cd ''' pwd ''';'];

%% 6 km all regions
KMRES=6;
RGNUMS=[1 2 3 4 5 6];
for avg_25hr = 0:1
    for region_num = RGNUMS
        matscript = ['master1(' num2str(region_num) ',' num2str(KMRES) ',' num2str(avg_25hr) ');" &'];
        matcmd = [matrun matvars matscript];
        disp(matcmd)
        eval(matcmd)
    end
end

%% 2 km no region 3
KMRES=2;
RGNUMS=[1 2 4 5 6];
for avg_25hr = 0:1
    for region_num = RGNUMS
        matscript = ['master1(' num2str(region_num) ',' num2str(KMRES) ',' num2str(avg_25hr) ');" &'];
        matcmd = [matrun matvars matscript];
        disp(matcmd)
        eval(matcmd)
    end
end

%% 1 km no region 3 or 6
KMRES=1;
RGNUMS=[1 2 4 5];
for avg_25hr = 0:1
    for region_num = RGNUMS
        matscript = ['master1(' num2str(region_num) ',' num2str(KMRES) ',' num2str(avg_25hr) ');" &'];
        matcmd = [matrun matvars matscript];
        disp(matcmd)
        eval(matcmd)
    end
end




%% To run a single region, resolution and 25hr avg or not
%
region_num              = 1;
radar_km_resolution     = 6;
avg_25hr                = 1;

master1(region_num,radar_km_resolution,avg_25hr)










