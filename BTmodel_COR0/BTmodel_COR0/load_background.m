

function [latitude,longitude,VORTclim,Fcor] = load_background(nlon,nlat,X,Y,MONname,model)
%  Usage [Psi,Phi] = init_case6(nlon,nlat,X,Y);
%  Input:  alpha  - rotation angle for test case 1
%          nlon   - number of longitudes
%          nlat   - number of latitudes
%          X      - meshgrid output (or slt_grid) of longitudes (lambda)
%          Y      - meshgrid output (or slt_grid) of latitudes (theta=asin(mu))
%  Output: 
%          Ucos   - u*cos velocity
%          Vcos   - v*cos velocity
%          Xi     - vorticity
%          Psi	  - stream function
%          Phi    - geopotential height (on unit sphere of radius a)  (nlon x nlat)
%	   Fcor   - Coriolis term in rotated coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon = X';  % note transpose so that lon corresponds to row dimension
lat = Y';

Omega   = 7.292e-5;  % units 1/s (earth angular velocity)

if strcmp(model,'ghg') || strcmp(model,'sst')
    model = 'control';
end

VORTclim=ncread(['inputs/vort200_' model '.nc'],'vorticity');
longitude=ncread(['inputs/vort200_' model '.nc'],'lon');
latitude=ncread(['inputs/vort200_' model '.nc'],'lat');

if strcmp(MONname,'JJA')
    VORTclim=squeeze(mean(VORTclim(:,:,6:8),3));
elseif strcmp(MONname,'M')
    VORTclim=VORTclim(:,:,5); 
elseif strcmp(MONname,'J')
    VORTclim=VORTclim(:,:,6); 
elseif strcmp(MONname,'J2')
    VORTclim=VORTclim(:,:,7); 
elseif strcmp(MONname,'A')
    VORTclim=VORTclim(:,:,8);  
elseif strcmp(MONname,'MJ')
    VORTclim=squeeze(mean(VORTclim(:,:,5:6),3));
elseif strcmp(MONname,'JJ')
    VORTclim=squeeze(mean(VORTclim(:,:,6:7),3));
elseif strcmp(MONname,'MJJ')
    VORTclim=squeeze(mean(VORTclim(:,:,5:7),3)); 
end

    
VORTclim=fliplr(VORTclim);
latitude=flipud(latitude);
Fcor = zeros(nlon,nlat);

for i=1:numel(lat)
  Fcor(i) = 2.0*Omega*sin(lat(i));
end
  
end
