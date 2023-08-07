% Douglas Cahl 2023
% University of South Carolina
%
% region_num = 1    region = 'USEGC'    % 1,2,6 km                              fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'
% region_num = 2    region = 'GAK'      % 1,2,6 km                              fullregionname = 'HFRADAR_US_Gulf_of_Alaska'
% region_num = 3    region = 'AKNS'     % 6 km                                  fullregionname = 'HFRADAR_Alaska_-_North_Slope'
% region_num = 4    region = 'USWC'     % 500m 1,2,6 km (500m is not supported)    fullregionname = 'HFRADAR_US_West_Coast'
% region_num = 5    region = 'USHI'     % 1,2,6 km                              fullregionname = 'HFRADAR_US_Hawaii'
% region_num = 6    region = 'PRVI'     % 2,6 km             fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands'
%
%% function master1_all.m
function [] = master1_all()

Nback       = 2; % files to go back when no data exists
sleeptime   = 0.1; % sleep between analysis in minutes
while true
    for radar_km_resolution = [1 2 6] % radar resolution (km)
        switch radar_km_resolution
            case 1
%                 RGNUMS = 1;
                RGNUMS=[1 2 4 5];
            case 2
%                 RGNUMS = 1;
                RGNUMS=[1 2 4 5 6];
            case 6
%                 RGNUMS = 1;
                RGNUMS= [1 2 3 4 5 6];
        end
        for avg_25hr = [0 1]
            for region_num = RGNUMS
                %%% region selection
                switch region_num
                    case 1
                        region = 'USEGC';    % 1,2,6 km
                        fullregionname = 'HFRADAR_US_East_and_Gulf_Coast';
                    case 2
                        region = 'GAK';      % 1,2,6 km
                        fullregionname = 'HFRADAR_US_Gulf_of_Alaska';
                    case 3
                        region = 'AKNS';     % 6 km
                        fullregionname = 'HFRADAR_Alaska_-_North_Slope';
                    case 4
                        region = 'USWC';     % 500m 1,2,6 km (500m is not supported)
                        fullregionname = 'HFRADAR_US_West_Coast';
                    case 5
                        region = 'USHI';     % 1,2,6 km
                        fullregionname = 'HFRADAR_US_Hawaii';
                    case 6
                        region = 'PRVI';     % 2,6 km
                        fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands';
                end

                % 25 hour average or hourly data
                switch avg_25hr
                    case 1
                        hr_name = '25hr';
                        hr_name_full = '25_Hour_Average';
                    otherwise
                        hr_name = 'hourly';
                        hr_name_full = 'Hourly';
                end


                % create data directory if missing
                dbase = 'data/'; % fiename in hours since 2019-09-29 20:00:00.000 UTC
                if ~exist(dbase,'dir')
                    mkdir(dbase) % create directory
                end

                % directory where hf currents are saved
                dout_data = [dbase 'analysis_data/']; % fiename in hours since 2019-09-29 20:00:00.000 UTC
                if ~exist(dout_data,'dir')
                    mkdir(dout_data) % create directory
                end
                dout_data = [dout_data  region  '/'];
                if ~exist(dout_data,'dir')
                    mkdir(dout_data) % create directory
                end
                dout_data = [dout_data hr_name '/'];
                if ~exist(dout_data,'dir')
                    mkdir(dout_data) % create directory
                end
                dout_data = [dout_data 'km' num2str(radar_km_resolution) '/'];
                if ~exist(dout_data,'dir')
                    mkdir(dout_data) % create directory
                end

                % analysis figures
                % output directory for saving analysis figures
                dout_fig = [dbase 'analysis_figs/'];
                if ~exist(dout_fig,'dir')
                    mkdir(dout_fig) % create directory
                end
                dout_fig = [dout_fig region '/'];
                if ~exist(dout_fig,'dir')
                    mkdir(dout_fig) % create directory
                end
                dout_fig = [dout_fig hr_name '/'];
                if ~exist(dout_fig,'dir')
                    mkdir(dout_fig) % create directory
                end
                dout_fig = [dout_fig 'km' num2str(radar_km_resolution) '/'];
                if ~exist(dout_fig,'dir')
                    mkdir(dout_fig) % create directory
                end

                % results directory for ftp upload
                tracks_dir = [dbase 'eddy_tracks/'];
                if ~exist(tracks_dir,'dir')
                    mkdir(tracks_dir) % create directory
                end
                tracks_dir = [tracks_dir region '/'];
                if ~exist(tracks_dir,'dir')
                    mkdir(tracks_dir) % create directory
                end
                tracks_dir = [tracks_dir hr_name '/'];
                if ~exist(tracks_dir,'dir')
                    mkdir(tracks_dir) % create directory
                end
                tracks_dir = [tracks_dir 'km' num2str(radar_km_resolution) '/'];
                if ~exist(tracks_dir,'dir')
                    mkdir(tracks_dir) % create directory
                end

                %    file1 = region '_eddy_tracks_' num2str(radar_km_resolution) 'km.csv'
                %    file2 = region '_eddy_tracks_' num2str(radar_km_resolution) 'km.npz'
                %    file3 = region '_eddy_tracks_' num2str(radar_km_resolution) 'km.mat'






                % while true

                % done in bash script now
                %    % create log directory if missing
                %    dlog = 'logs/' % fiename in hours since 2019-09-29 20:00:00.000 UTC
                %    logfn = region num2str(radar_km_resolution) hrname '.txt'
                %    if not path.isdir(dlog):
                %        mkdir(dlog) % create directory

                % url index
                url1 = ['http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/' region ...
                    '/' num2str(radar_km_resolution) ...
                    'km/' hr_name ...
                    '/RTV/' fullregionname ...
                    '_' num2str(radar_km_resolution) ...
                    'km_Resolution_' hr_name_full ...
                    '_RTV_best.ncd.html'];

                % get data set lat lon size
                try
                    mystr = urlread(url1);
                catch
                    disp('website not giving description')
                    continue
                end
                if isempty(mystr)
                    continue
                end

                % time offset
                % disp(mystr)
                i_str = strfind(mystr,'units: hours since ');
                if isempty(i_str)
                    continue
                end
                mystr1_time = mystr(i_str(1)+19:i_str(1)+45);
                %    disp(mystr1_time) % offset from date
                % hours
                i_str = strfind(mystr,'Float64 time_run[time =');
                mystr1 = mystr(i_str+24:i_str+31);
                i_str = strfind(mystr1,'];');
                tnum = str2double(mystr1(1:i_str-1)); % latest time + 1, to set numfiles for this loop
                %    nc_time = tnum - 1;
                yyyy = str2double(mystr1_time(1:4));
                mm = str2double(mystr1_time(6:7));
                dd = str2double(mystr1_time(9:10));
                hh = str2double(mystr1_time(12:13));

                now = datenum(datetime('now','timezone','utc'));
                start = datenum(yyyy,mm,dd,hh,0,0);
                hours_now = (now - start)*24; % change to hours_now = hours_now - 20-100;

                % to setup lat and nom
                i_str = strfind(mystr,'Float32 lat[lat = ');
                mystr1 = mystr(i_str(1)+18:i_str(1)+29);
                i_str = strfind(mystr1,'];');
                latnum = str2double(mystr1(1:i_str-1));

                i_str = strfind(mystr,'Float32 lon[lon = ');
                mystr1 = mystr(i_str(1)+18:i_str(1)+29);
                i_str = strfind(mystr1,'];');
                lonnum = str2double(mystr1(1:i_str-1));



                %    now = datetime.datetime.utcnow()
                %    start = datetime.datetime(2019,9,29)
                %    diff = now - start
                %    days, seconds = diff.days, diff.seconds
                %    hours_now = days * 24 + seconds // 3600
                %    hours_now = hours_now - 20; % change to hours_now = hours_now - 20-100;



                disp(['current time ahead by ' num2str(round(hours_now,1)) ' hours'])
                if avg_25hr == 1
                    disp('25 hr average 12.5 hours behind')
                end

                files = dir(dout_data);
                if length(files) < 3
                    numfiles = Nback;
                    nc_time_start = tnum - Nback;
                    disp(['last data file: ' num2str(tnum)])
                    disp(['first run, starting ' num2str(Nback) ' hours back from last data'])
                else
                    lastfile = files(end).name;
                    lasthr = strsplit(lastfile,'_');
                    lasthr = lasthr{end};
                    lasthr = strsplit(lasthr,'.');
                    lasthr = str2double(lasthr{1});
                    nc_time_start = lasthr + 1;
                    numfiles = floor(hours_now - lasthr);
                    disp(['last file +1 = ' num2str(nc_time_start)])
                    disp(['last data file: ' num2str(tnum)])
                    disp(['catchup num  = ' num2str(numfiles)])
                end

                %    numfiles = 1


                %     return
                eddy1_identification;
                eddy2_track_eddies()
                eddy3_analyze_eddy_tracks()

                %    break

                % uploads results to server
                % disp('sending results')
                % try:
                %     session = ftplib.FTP('ftp_ip','ftpuser','ftp_pass') % ftp server
                %     file = open(fdir file1,'rb')                  % file to send
                %     session.storbinary('STOR ' file1, file)       % send the file
                %     file = open(fdir file2,'rb')
                %     session.storbinary('STOR ' file2, file)
                %     file = open(fdir file3,'rb')
                %     session.storbinary('STOR ' file3, file)
                %     file.close()                                    % close file and FTPsession = ftplib.FTP('104.128.235.148','ftp_home','homehome')
                %     session.quit()
                % except:
                %     disp('ftp server down, no results going to webserver')


                disp(['sleep for ' num2str(sleeptime) ' minutes ...'])
                pause(60*sleeptime) % sleep until next run
                %  end


            end
        end
    end
end

