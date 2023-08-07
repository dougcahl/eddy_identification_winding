# eddy_identification_winding
Eddy identificationfor HF radar data sets using the Winding Angle method using matlab
Modified and written by Douglas Cahl, PhD 
University of South Carolina 2023

Add m_map and the etopo1 data inside m_map/etopo1/etopo1_ice_g_i2.bin, etopo1_ice_g_i2.hdr

Use <code>eddy_uvdata.m</code> for your own lon,lat,u,v data set or for single netcdf files

Use <code>eddy_uvdata_loop.m</code> for your own lon,lat,u,v, time data set

Use <code>eddy_tracking.m</code> tracks the eddies found from above script

Use <code>analyze_eddy_tracks.m</code> plots the eddy tracks found from above script and some data analysis


Use <code>eddy_netcdf.m</code> to loop through large downloaded data sets of netcdf files in 1,2 or 6km resolution

More information available at http://douglascahl.com/eddy/

# all regions
Eddy identification for HF radar data sets using the Winding Angle method. Supports all 6 US regions of HFR data for downloading and processing in 1,2 and 6 km resolution. You can run all the regions and resolutions simultaneously with enough ram. 


see <code>master1_all.m</code> for details on how to run


<p>
</p>
A video description, data download and an interactive realtime eddy analysis overlayed on google maps is available at http://douglascahl.com/eddy/
