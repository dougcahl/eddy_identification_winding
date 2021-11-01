#!/bin/bash

# cuts header line of all csv files
#find . | grep 'km.csv' | xargs -I {} sed 1d "{}" > "{}1"
#find . | grep 'avg25hr_1_eddy_tracks_6km.csv' | xargs -I {} sed 1d "{}" > "{}1"

BD="/home/doug/eddy_identification_winding/Data/eddy_tracks/"

find ${BD} | grep 'avg25hr_1_eddy_tracks_6km.csv' | tail -1 | xargs -I{} head -1 {} > ${BD}eddy_tracks_avg25hr_1_6km.csv
find ${BD} | grep 'avg25hr_0_eddy_tracks_6km.csv' | tail -1 | xargs -I{} head -1 {} > ${BD}eddy_tracks_avg25hr_0_6km.csv
find ${BD} | grep 'avg25hr_1_eddy_tracks_2km.csv' | tail -1 | xargs -I{} head -1 {} > ${BD}eddy_tracks_avg25hr_1_2km.csv
find ${BD} | grep 'avg25hr_0_eddy_tracks_2km.csv' | tail -1 | xargs -I{} head -1 {} > ${BD}eddy_tracks_avg25hr_0_2km.csv
find ${BD} | grep 'avg25hr_1_eddy_tracks_1km.csv' | tail -1 | xargs -I{} head -1 {} > ${BD}eddy_tracks_avg25hr_1_1km.csv
find ${BD} | grep 'avg25hr_0_eddy_tracks_1km.csv' | tail -1 | xargs -I{} head -1 {} > ${BD}eddy_tracks_avg25hr_0_1km.csv

find ${BD} | grep 'avg25hr_1_eddy_tracks_6km.csv' | xargs -I{} sed 1d {} >> ${BD}eddy_tracks_avg25hr_1_6km.csv
find ${BD} | grep 'avg25hr_0_eddy_tracks_6km.csv' | xargs -I{} sed 1d {} >> ${BD}eddy_tracks_avg25hr_0_6km.csv
find ${BD} | grep 'avg25hr_1_eddy_tracks_2km.csv' | xargs -I{} sed 1d {} >> ${BD}eddy_tracks_avg25hr_1_2km.csv
find ${BD} | grep 'avg25hr_0_eddy_tracks_2km.csv' | xargs -I{} sed 1d {} >> ${BD}eddy_tracks_avg25hr_0_2km.csv
find ${BD} | grep 'avg25hr_1_eddy_tracks_1km.csv' | xargs -I{} sed 1d {} >> ${BD}eddy_tracks_avg25hr_1_1km.csv
find ${BD} | grep 'avg25hr_0_eddy_tracks_1km.csv' | xargs -I{} sed 1d {} >> ${BD}eddy_tracks_avg25hr_0_1km.csv



#find . | grep 'km.csv' | xargs -I{} sed 1d {} > eddy_tracks.csv
