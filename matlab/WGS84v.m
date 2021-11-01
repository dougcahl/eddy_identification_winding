%% WGSV84v.m
function [d,phi] = WGS84v(lat,lon,lat_o,lon_o)
% [d,phi] = WGSD84v(lat,lon,lat_o,lon_o)
%
% Function to compute distance and angles between points on the WGS-84 
% ellipsoidal Earth to within a few millimeters of accuracy  using Vincenty's algorithm.
%
%  T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
%  with Application of Nested Equations", Survey Review, vol. 23, no. 176,
%  April 1975, pp 88-93
%
%% Input:
%  lat   = Geodetic latitude of grid point(s) (degrees) 
%  lon   = Geodetic longitude of grid point(s) (degrees)
%  lat_o = Geodetic latitude of Radar site (degrees) 
%  lon_o = Geodetic longitude of Radar site (degrees)
%
%  (lat, lon) can be vectors or matrices
%  (lat_o,lon_o) must be a single point (radar site)
%
%% Output:
%  d   = distance between each point (lat,lon) and (lat_o,lon_o) (in km)
%  phi = azimuthal angle (in degs) from radar site to grid point(s)
%               
%% Example:
%  [d,phi]=WGS84v( [25.2804153,25.4604010],[-80.1641749,-79.0696098],25.713333,-80.150833)
%
%  s =    47.9775 km
%        112.1807 km
%  phi = 181.6050 degs North
%        104.2304 degs North
%
%% Credits
%  Based on an original code by Michael Kleder
%  downloaded from: 
%  http://www.mathworks.com/matlabcentral/fileexchange/5379
%  
%  Modified for Radar Application and Vectorized by George Voulgaris
%  Department of Earth & Ocean Sciences
%  University of South Carolina, 
%  Columbia, SC 29208, USA
%
%  Notes:(1) Error correcting code, convergence failure traps, antipodal corrections,
%            polar error corrections, WGS84 ellipsoid parameters, testing, and comments
%            written by Michael Kleder, 2004.
%        (2) Vincenty describes his original algorithm as precise to within
%            0.01 millimeters, subject to the ellipsoidal model.
%        (3) Essentially antipodal points are treated as exactly antipodal,
%            potentially reducing accuracy by a small amount.
%        (4) Failures for points exactly at the poles are eliminated by
%            moving the points by 0.6 millimeters.
%        (5) Vincenty's azimuth formulas are not implemented in this
%            version, but are available as comments in the code.
%        (6) The Vincenty procedure was transcribed verbatim by Peter Cederholm,
%            August 12, 2003. It was modified and translated to English by Michael Kleder.
%            Mr. Cederholm's website is http://www.plan.aau.dk/~pce/
%        (7) Code to test the disagreement between this algorithm and the
%            Mapping Toolbox spherical earth distance function is provided
%            as comments in the code. The maximum differences are:
%            Max absolute difference: 38 kilometers
%            Max fractional difference: 0.56 percent
%
%% Copyright 2019, George Voulgaris, University of South Carolina
%
% This file is part of matWERA.
%
% matWERA is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% If you find an error please notify G. Voulgaris at gvoulgaris@geol.sc.edu
%
%%
[m,n] = size(lon);     % Check dimensions of input data
lon1 = lon(:);
lat1 = lat(:);
%
% -----------------------Input check --------------------------------------
%
ier = find((abs(lat1)>90 | abs(lat_o)>90),1);
if ~isempty(ier)
    error('Some input latitudes where found to be > +/- 90 degrees.')
end
%  WGS84 earth ellipsoid axis lengths in meters:
a = 6378137;          % definitionally
b = 6356752.31424518; % computed from WGS84 earth flattening coefficient definition
% convert inputs in degrees to radians:
lat1 = lat1  * 0.0174532925199433;
lon1 = lon1  * 0.0174532925199433;
lat2 = lat_o * 0.0174532925199433;
lon2 = lon_o * 0.0174532925199433;
% -- correct for errors at exact poles by adjusting 0.6 millimeters -------
i1 = find(abs(pi/2-abs(lat1)) < 1e-10);
if ~isempty(i1)
    lat1(i1) = sign(lat1(i1))*(pi/2-(1e-10));
end
i2 = find(abs(pi/2-abs(lat2)) < 1e-10);
if ~isempty(i2)
    lat2(i2) = sign(lat2(i2))*(pi/2-(1e-10));
end
% -------------------------------------------------------------------------
f    = (a-b)/a;
U1   = atan((1-f)*tan(lat1));
U2   = atan((1-f)*tan(lat2));
lon1 = mod(lon1,2*pi);
lon2 = mod(lon2,2*pi);
L    = abs(lon2-lon1);
i3   = find(L > pi);
if ~isempty(i3)
    L(i3) = 2*pi - L(i3);
end
lambda    = L;
lambdaold = zeros(size(lambda));
itercount = 0; 

%while ~itercount | abs(lambda-lambdaold) > 1e-12  % force at least one execution
while itercount<=50 && max(abs(lambda-lambdaold)) > 1e-12
    itercount = itercount+1;
%     if itercount > 50
%          warning('Points are essentially antipodal. Precision may be reduced slightly.');
%         lambda = pi;
%         break
%     end
    lambdaold  = lambda;
    sinsigma   = sqrt((cos(U2).*sin(lambda)).^2+(cos(U1).*...
                 sin(U2)-sin(U1).*cos(U2).*cos(lambda)).^2);
    cossigma   = sin(U1).*sin(U2)+cos(U1).*cos(U2).*cos(lambda);
    sigma      = atan2(sinsigma,cossigma);
    alpha      = asin(cos(U1).*cos(U2).*sin(lambda)./sin(sigma));
    cos2sigmam = cos(sigma)-2*sin(U1).*sin(U2)./cos(alpha).^2;
    C          = f/16*cos(alpha).^2.*(4+f*(4-3*cos(alpha).^2));
    lambda     = L+(1-C)*f.*sin(alpha).*(sigma+C.*sin(sigma).*...
                (cos2sigmam+C.*cos(sigma).*(-1+2*cos2sigmam.^2)));
end
%  correct for convergence failure in the case of essentially antipodal points
io = find (lambda>pi);    % In case of antipodal points
if ~isempty(io)
    lambda(io)=lambda(io)*0+pi;
end
u2         = cos(alpha).^2*(a^2-b^2)/b^2;
A          = 1+u2/16384.*(4096+u2.*(-768+u2.*(320-175*u2)));
B          = u2./1024.*(256+u2.*(-128+u2.*(74-47*u2)));
deltasigma = B.*sin(sigma).*(cos2sigmam+B/4.*(cos(sigma).*(-1+2*cos2sigmam.^2)...
            -B/6.*cos2sigmam.*(-3+4*sin(sigma).^2).*(-3+4*cos2sigmam.^2)));
s          = b*A.*(sigma-deltasigma);
%
%==========================================================================
% Vicenty's azimuth calculation code (results in radians)
% ----------- From point #1 (Grid) to point #2 (Radar)---------------------
%a12     = atan2(cos(U2).*sin(lambda),cos(U1).*sin(U2)-sin(U1).*cos(U2).*cos(lambda));
%i3      = find(a12 < 0);
%a12(i3) = a12(i3)+pi; 
% ----------- from point #2 (Radar) to point #1 (Grid)---------------------
a21     = atan2(cos(U1).*sin(lambda),-sin(U1).*cos(U2)+cos(U1).*sin(U2).*cos(lambda));
i4      = find(a21 < 0);
a21(i4) = a21(i4)+pi; 
% ----------- Cases when Longitude difference (L) is 0<L<pi--------------------
i5      = find((L>0) & (L<pi));
a21(i5) = a21(i5) + pi;
i6      = find((lon2-lon1)<0);
%-Make sure that the azimuth angle is always on the direction from radar to grid point -
a21(i6) = 2*pi - a21(i6);     
%--------------------------------------------------------------------------
a21     = a21*180/pi;       % Convert to degrees
s       = s/1000;           % Convert to km
phi     = reshape(a21,m,n); % Convert back to original size
d       = reshape(s,m,n);
end