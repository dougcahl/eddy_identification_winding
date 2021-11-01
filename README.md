# eddy_identification_winding
Eddy identificationfor HF radar data sets using the Winding Angle method

<code>
call with master1.py avg_25hr radar_km_resolution region_num
where 
region_num = 1    region = 'USEGC'    # 1,2,6 km                              fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'
region_num = 2    region = 'GAK'      # 1,2,6 km                              fullregionname = 'HFRADAR_US_Gulf_of_Alaska'
region_num = 3    region = 'AKNS'     # 6 km                                  fullregionname = 'HFRADAR_Alaska_-_North_Slope'
region_num = 4    region = 'USWC'     # 500m 1,2,6 km (500m is not supported)    fullregionname = 'HFRADAR_US_West_Coast'
region_num = 5    region = 'USHI'     # 1,2,6 km                              fullregionname = 'HFRADAR_US_Hawaii'
region_num = 6    region = 'PRVI'     # 2,6 km             fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands'
example: pythong3 master1.py 1 6 1

run with:   ./eddy_all_regions.sh > master.log
find with:  sudo ps -af | grep master1
kill with:  sudo pkill -f 'master1.py'
check logs: ls *.txt | xargs -I{} tail -1 {}
</code>

A video description, data download and an interactive realtime eddy analysis overlayed on google maps is available at http://douglascahl.com/eddy/
