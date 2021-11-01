#!/bin/bash

mkdir -p logs

# call with master1.py avg_25hr radar_km_resolution region_num
# where 
# region_num = 1    region = 'USEGC'    # 1,2,6 km                              fullregionname = 'HFRADAR_US_East_and_Gulf_Coast'
# region_num = 2    region = 'GAK'      # 1,2,6 km                              fullregionname = 'HFRADAR_US_Gulf_of_Alaska'
# region_num = 3    region = 'AKNS'     # 6 km                                  fullregionname = 'HFRADAR_Alaska_-_North_Slope'
# region_num = 4    region = 'USWC'     # 500m 1,2,6 km (500m is not supported)    fullregionname = 'HFRADAR_US_West_Coast'
# region_num = 5    region = 'USHI'     # 1,2,6 km                              fullregionname = 'HFRADAR_US_Hawaii'
# region_num = 6    region = 'PRVI'     # 2,6 km             fullregionname = 'HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands'
# example: pythong3 master1.py 1 6 1

# run with:   ./eddy_all_regions.sh > master.log
# find with:  sudo ps -af | grep master1
# kill with:  sudo pkill -f 'master1.py'
# check logs: ls *.txt | xargs -I{} tail -1 {}


# 6 km all regions
KMRES=6
RGNUMS=(1 2 3 4 5 6)
for avg_25hr in {0..1}
do
    for region_num in "${RGNUMS[@]}"
    do
#        cmd="python3 -u master1.py ${avg_25hr} ${KMRES} ${region_num} >> logs/avg25hr_${avg_25hr}_res_${KMRES}_region_${region_num}.txt 2>&1 &"
        cmd="python3 -u master1.py ${avg_25hr} ${KMRES} ${region_num} >> logs/avg25hr_${avg_25hr}_res_${KMRES}_region_${region_num}.txt &"
        echo $cmd
        eval $cmd
    done
done

# # 2 km no region 3
# KMRES=2
# RGNUMS=(1 2 4 5 6)
# for avg_25hr in {0..1}
# do
#     for region_num in "${RGNUMS[@]}"
#     do
# #        cmd="python3 -u master1.py ${avg_25hr} ${KMRES} ${region_num} >> logs/avg25hr_${avg_25hr}_res_${KMRES}_region_${region_num}.txt 2>&1 &"
#         cmd="python3 -u master1.py ${avg_25hr} ${KMRES} ${region_num} >> logs/avg25hr_${avg_25hr}_res_${KMRES}_region_${region_num}.txt &"
#         echo $cmd
#         eval $cmd
#     done
# done

# # 1 km no region 3 or 6
# KMRES=1
# RGNUMS=(1 2 4 5)
# for avg_25hr in {0..1}
# do
#     for region_num in "${RGNUMS[@]}"
#     do
# #        cmd="python3 -u master1.py ${avg_25hr} ${KMRES} ${region_num} >> logs/avg25hr_${avg_25hr}_res_${KMRES}_region_${region_num}.txt 2>&1 &"
#         cmd="python3 -u master1.py ${avg_25hr} ${KMRES} ${region_num} >> logs/avg25hr_${avg_25hr}_res_${KMRES}_region_${region_num}.txt &"
#         echo $cmd
#         eval $cmd
#     done
# done



#
#
#for avg_25hr in range(1,7):
#    for radar_km_resolution in range(1,7):
#        for region_num in range(1,7):
#
#
#
#

#
# To run a single region and resolution
#
#python3 -u master1.py 0 6 1 >> logs/eddy_master_usegc6_hourly.txt &
#python3 -u master1.py 0 2 1 >> logs/eddy_master_usegc2_hourly.txt &
#python3 -u master1.py 0 1 1 >> logs/eddy_master_usegc1_hourly.txt &
#
#python3 -u master1.py 1 6 1 >> logs/eddy_master_usegc6_25hr.txt &
#python3 -u master1.py 1 2 1 >> logs/eddy_master_usegc3_25hr.txt &
#python3 -u master1.py 1 1 1 >> logs/eddy_master_usegc1_25hr.txt &


#python3 -u eddy_master_usegc6_hourly.py >> logs/eddy_master_usegc6_hourly.txt &
#python3 -u eddy_master_usegc2_hourly.py >> logs/eddy_master_usegc2_hourly.txt &
#python3 -u eddy_master_usegc1_hourly.py >> logs/eddy_master_usegc1_hourly.txt &
#
#python3 -u eddy_master_usegc6_25hr.py >> logs/eddy_master_usegc6_25hr.txt &
#python3 -u eddy_master_usegc2_25hr.py >> logs/eddy_master_usegc2_25hr.txt &
#python3 -u eddy_master_usegc1_25hr.py >> logs/eddy_master_usegc1_25hr.txt &

#pid1 = subprocess.Popen([sys.executable, "eddy_master_usegc6_hourly.py"])
#pid2 = subprocess.Popen([sys.executable, "eddy_master_usegc2_hourly.py"])
#pid3 = subprocess.Popen([sys.executable, "eddy_master_usegc1_hourly.py"])
#pid4 = subprocess.Popen([sys.executable, "eddy_master_usegc6_25hr.py"])  
#pid5 = subprocess.Popen([sys.executable, "eddy_master_usegc2_25hr.py"])
#pid6 = subprocess.Popen([sys.executable, "eddy_master_usegc1_25hr.py"])

