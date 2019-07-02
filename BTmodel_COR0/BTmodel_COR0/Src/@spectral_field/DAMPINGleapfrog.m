
function Vstep = DAMPINGleapfrog(VORTstep_t,DIFForder,DIFFcoeff,dt) 
% DEL2 - Laplacian of a spectral_grid field.
% Usage: g = del2(f)
% Compute the Laplacian from the grid point values using the spectral
% transform and modifying the spectral coefficients.  The grid point
% values are returned from a spherical harmonic transform synthesis
%--------------------------------------------------------
%  Input: 
%    f - spectral_field class objects
%        In particular:
%           f.G  - field gauss_grid
%           f.gp - field grid point values
%           f.sc - field spectral coefficients
%  Output:
%    g - spectral_field object = del2(f)
%  Local
%    xf - matrix of Fourier coefficients ordered (m,j)
%--------------------------------------------------------

mm = get(VORTstep_t.G,'mm');
nn = get(VORTstep_t.G,'nn');
kk = get(VORTstep_t.G,'kk');
aradius= get(VORTstep_t.G,'radius');
%--------------------------------------------------------
% transform the fields to spectral space
Vstep = shtrana(VORTstep_t);    

Vstep_sc = get(Vstep,'sc');

% modify the spectral coefficients
for n = 1:nn
  	m = 0:n;
    SIGMA=((n+1)*(n))/(aradius^2);
    Denom=(-1)*DIFFcoeff*((-1)^DIFForder)*(2*dt)*(SIGMA^DIFForder) + 1;
	Vstep_sc(n+1,m+1) = Vstep_sc(n+1,m+1)/Denom;  %is correct
end
Vstep_sc(1,1) = complex(0.0,0.0);
Vstep = set(Vstep,'sc',Vstep_sc);% set the new coefficients
Vstep = shtrans(Vstep);    % inverse transform to get new grid point values 
% all done  - f is unchanged, but new spectral_field g created.
