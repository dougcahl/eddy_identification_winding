function [X,Y,UTMzone]=geog2utm_nodisp(lon,lat,LON,LAT)
%
%  Convert Lat and Lon to UTM (Elipsoid) system (in Km) with Radar site being (0,0). 
%  
%    function [X,Y]=geog2utm(longitude,latitude)
%
%  Function to convert geographic coordinates (Longitude &
%  latitude) to coordinates on the Universal Transverse
%  Mercator (Elipsoid) system. The program finds automatically 
%  the UTM zone from the center of latitudes.
%
%  function [X,Y]=geog2utm(lon,lat,LON,LAT)
%
%  as above but instead of using the UTM zones the user defines 
%  the origin coordinates LAT and LON.
%
%  INPUT : 	lon, lat : longitude and latitude in degrees and decimals
%                          of a degree (South & West must be negative, 
%                          i.e. 5.2N 3.3W is +5.2 -3.3).
%		LON, LAT : [OPTIONAL] Latitude and longitude of the origin poin 
%                          (0.0,0.0)(South & West must be negative).
%  OUTPUT:	X  , Y   : Eastings and Northings (in km) from the origin point 
%                          
%  J.P. Snyder, 1987, Map Projections - A Working Manual. 
%  USGS Professional Paper 1395.
%
%  G Voulgaris, USC  FEB. 2000
% 
%  Revision 1, Jan. 2011: Took care of the UTM 500km offset automatically.
%
echo off
%
lat     = lat(:);
lon     = lon(:);
UTMzone = 0;
%
if (nargin > 4 || nargin==1 || nargin==3)
   error('geog2utm.m: Only 2 or 4 arguments are allowed')
end
% Find the UTM Zone and corresponding Meridian using the mean longitude
if nargin == 2 
    LL          = min(lon)+(max(lon)-min(lon))/2;
    UTMrange    = (-180:6:180)';
    UTMmeridian = (-177:6:178)';
    for i = 1:60;
      if (UTMrange(i)<LL && UTMrange(i+1)>=LL); 
          UTMzone=i; 
      end
    end
    XOFF = 500.0;   % The 500Km offset required for UTM coordinates
    LAT  = 0;
    LON  = UTMmeridian(UTMzone);
else
    XOFF = 0.0;       % In case of local coordinate system with given origin
end
%
% Constants
%
% disp([' > UTM Zone:',num2str(UTMzone),' ORIGIN: ',num2str(LON),'E, ',num2str(LAT),'N <'])
%
alpha   = 6378206.4;
e2      = 0.00676866;
DEG2RAD = (2*pi/360);
ko      = 0.9996;
%
% Converting to rads
%
lat = lat*DEG2RAD;
lon = lon*DEG2RAD;
LAT = LAT*DEG2RAD;
LON = LON*DEG2RAD;
%
ee2 = e2/(1-e2);
N   = alpha./sqrt((1-e2*sin(lat).^2));
T   = tan(lat).^2;
C   = ee2*cos(lat).^2;
A   = (lon-LON).*cos(lat);
M   = 111132.0894*lat/DEG2RAD-16216.94*sin(2*lat)+17.21*sin(4*lat)-0.02*sin(6*lat);
Mo  = 111132.0894*LAT/DEG2RAD-16216.94*sin(2*LAT)+17.21*sin(4*LAT)-0.02*sin(6*LAT);
%
X   = ko*N.*(A+(1-T+C).*(A.^3/6)+(5-18*T+T.^2+72*C-58*ee2).*(A.^5/120));
Y   = ko*(M-Mo+N.*tan(lat).*((A.^2/2)+(5-T+9*C+4*C.^2).*(A.^4/24)+(61-58*T+T.^2+600*C-330*ee2).*(A.^6/720)));
%
X   = X/1000+XOFF; % in Km
Y   = Y/1000;      % in Km
end