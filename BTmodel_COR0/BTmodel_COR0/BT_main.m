function BT_main(LONGoutput)%, LONraster)%, model, MONname)

GETinitTEND=0;
RUNstationary=1;
FORCINGon=1;
FORCEmon='M';
MONname='J'
model='rpmcontrol';

%MONname=strcat(MONname);
%model=strcat(model);


OUTfname = ['AAA' model '_' MONname '_clim_' FORCEmon '_forcing_NEW']
%OUTfname = [model '_' MONname '_clim_' int2str(LONraster)]

FORCEfname=['outputs/INITtend_' model '_' MONname 'clim' '.mat'];

if LONGoutput==1
    OUTfname = [OUTfname '_LONGout']
end 

PLOTfig=1;

FORCINGlinearFACTOR=1e-6;

NONDIVadjust=0;

WHITEnoise=0;
WHITEamp=FORCINGlinearFACTOR*(1e-11);

GAUSSforcing=1;

g1=1;
LATpeak=-10;
LONpeak=62.5;
LATsigma=5;
LONsigma=10;

g2=0;
LATpeak2=7.5;
LONpeak2=220;
LATsigma2=2.5;
LONsigma2=10;

g3=0;
LATpeak3=7.5;
LONpeak3=305;
LATsigma3=2.5;
LONsigma3=10;

GAUSScoeff=FORCINGlinearFACTOR*(2.5e-11);

COSforcing=0;
COSforcingsingle=0;

COScoeff=FORCINGlinearFACTOR*(2.5e-11);


RWSforcing=0;
RWSfile=['inputs/forcing_' model '.nc'];
SOURCErws=ncread(RWSfile,'forcing');

if strcmp(FORCEmon,'JJA')
    SOURCErws=squeeze(mean(SOURCErws(:,:,6:8),3));
elseif strcmp(FORCEmon,'M')
    SOURCErws=SOURCErws(:,:,5); 
elseif strcmp(FORCEmon,'J')
    SOURCErws=SOURCErws(:,:,6); 
elseif strcmp(FORCEmon,'J2')
    SOURCErws=SOURCErws(:,:,7); 
elseif strcmp(FORCEmon,'A')
    SOURCErws=SOURCErws(:,:,8);  
elseif strcmp(FORCEmon,'MJ')
    SOURCErws=squeeze(mean(SOURCErws(:,:,5:6),3));
elseif strcmp(FORCEmon,'JJ')
    SOURCErws=squeeze(mean(SOURCErws(:,:,6:7),3));
elseif strcmp(FORCEmon,'MJJ')
    SOURCErws=squeeze(mean(SOURCErws(:,:,5:7),3)); 
end

LATrws=ncread(RWSfile,'lat');
LONrws=ncread(RWSfile,'lon');

%==== SELECT SOURCE REGION ====
SOURCErws(:,LATrws>10)=0;
SOURCErws(:,LATrws<0)=0;
% IND + ATL
%SOURCErws(33:107,:)=0;
%SOURCErws(LONrws<45,:)=0;
%IND
%SOURCErws(LONrws<45,:)=0;
%SOURCErws(LONrws>90,:)=0;
%ATL
%SOURCErws(LONrws<300,:)=0;
%NIND + NATL
%SOURCErws(17:33,:)=0;
%SOURCErws(LONrws>300,:)=0;
% PAC
%SOURCErws(LONrws<220)=0;
%SOURCErws(LONrws>250)=0;

%SOURCErws(LONrws<270,:)=0;
%SOURCErws(LONrws>330,:)=0;
%SOURCErws(33:96,:)=0;
%SOURCErws(LONrws<90, LATrws>10)=0;
%SOURCErws(LONrws>300, LATrws>10)=0;
%SOURCErws(96:107,:)=0;


SOURCErws=flipdim(SOURCErws,2);




%DIFFUSIONon=0;

DIFForder=2;
DIFFcoeff=-5e16;%2.4 for JJA

DRAGcoeff=1.2e-6;
%DRAGcoeff=0;

%ROBERTcoeff=0.15;


if RUNstationary==1
load(FORCEfname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Begin method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
string = 'Begin'



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control initialization  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nj = input('enter the value of N   ')
%ndays = input('enter the number of days to integrate ')

%-- T42
nj=64;     % number of latitudes
ndays = 25;  % Number of days to integrate the model
T = 86400*ndays   % time interval of the solution in seconds
%nsteps=5;  % number of steps to take
%dt = T/nsteps;   % model time step
dt = 900   % model time step (15 minutes)
nsteps = T/dt
nout=96;    % frequency of model output in number of steps
if LONGoutput==1
    nout=4;
end
t = 0:dt:T;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U250=[];
V250=[];
psi250=[];
vort250=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grid initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain Spherical grid (ni,nj)
G = gauss_grid('T42',nj);
aradius = 6.37122e6; % units - m (earth radius)
G = set(G,'radius',aradius);  % set the Earths radius for the grid operators
ni = get(G,'ni');
xg = get(G,'xg');
yg = get(G,'yg');
%setup the Slt extended grid
A = slt_grid(xg,yg);
LON = get(A,'X');
LAT = get(A,'Y');
LAT = asin(LAT);  % adjust to longitude for initialization





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Initial Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define the spectral_fields on the spectral_grid G
%Prognostic:
%ETA= spectral_field('Potential Vorticity',G);
%Diagnostic
PSI_t= spectral_field('Streamfunction (t)',G);

VORT_tp1= spectral_field('Vorticity (t+1)',G);
VORT_t= spectral_field('Vorticity (t)',G);
VORT_tm1= spectral_field('Vorticity (t-1)',G);

DELTA= spectral_field('Divergence (t)',G);
Z= zeros(ni,nj);
FORCING = zeros(ni,nj);    
SOURCE = zeros(ni,nj);
SOURCEcos = zeros(ni,nj);   


Uadvect_t = spectral_field('Zonal vorticity flux (t)',G);
Vadvect_t = spectral_field('Meridional vorticity flux (t)',G);

VORTstep_t = spectral_field('Vorticity (t) - dt*ADVECT_t',G);

DIFFUSE_t = spectral_field('Vorticity (t) - dt*ADVECT_t',G);

%WHITEsmooth = spectral_field('Smoothed white noise forcing',G);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializiation: U250 JJA CLIMATOLOGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[lat,lon,PSI_gp_INIT,Fcor] = init_case_sf250_1979to2016(ni,nj,LON,LAT); % gp is PSI, gp1 is PHI
%[lat,lon,PSI_gp_INIT,Fcor] = init_case_sf500month_1979to2016(ni,nj,LON,LAT,MONname);

%[lat,lon,VORT_gp_INIT,Fcor] = init_case_vortLEVELmonth_1979to2016(ni,nj,LON,LAT,MONname,LEVELtype,LEVELnum);

[lat,lon,VORT_gp_INIT,Fcor] = load_background(ni,nj,LON,LAT,MONname,model);


VORT_t=set(VORT_t,'gp',VORT_gp_INIT);
%FCOR = set(VORT,'gp',Fcor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET GAUSSIAN FORCING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FORCINGon==1
    
    if WHITEnoise==1
        FORCING = WHITEamp*randn(size(FORCING));
        for jj=1:nj
            FORCING(:,jj) = FORCING(:,jj) - mean(FORCING(:,jj));
        end 
        
%         WHITEsmooth = set(WHITEsmooth,'gp',FORCING);
%         WHITEsmooth = shtrana(WHITEsmooth);
%         WHITEsmooth = del2inv(WHITEsmooth);
   end

    if GAUSSforcing==1
        for ii=1:128
            for jj=1:64
                SOURCE(ii,jj)=g1*exp(-1*(((lon(ii)-LONpeak)^2)/(LONsigma^2) + ((lat(jj)-LATpeak)^2)/(LATsigma^2)))+g2*exp(-1*(((lon(ii)-LONpeak2)^2)/(LONsigma2^2) + ((lat(jj)-LATpeak2)^2)/(LATsigma2^2)))+g3*exp(-1*(((lon(ii)-LONpeak3)^2)/(LONsigma3^2) + ((lat(jj)-LATpeak3)^2)/(LATsigma3^2)));
            end
        end
        %SOURCE=GAUSScoeff*SOURCE.*Fcor;
        SOURCE=GAUSScoeff*SOURCE;
        FORCING = FORCING + SOURCE;
    end
    
    if RWSforcing==1
        FORCING = FORCING + FORCINGlinearFACTOR*SOURCErws;
    end
    
    if COSforcing==1
        if mod(LATraster,2)==0
            LONraster = LONraster+.5;
        end
        if mod(LATraster,2) ~= 0 && LONraster == 1
            for ii=LONraster*4-3:LONraster*4-3+1
                for jj=LATraster+33-2:LATraster+33+2
                    SOURCEcos(ii,jj)=((cos(0.5*pi*(lon(ii)-lon(LONraster*4-3))/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(LATraster+33))/(lat(31)-lat(33))))^2);
                end
            end
            for ii=127:128
                for jj=LATraster+33-2:LATraster+33+2
                    SOURCEcos(ii,jj)=SOURCEcos(ii,jj)+((cos(0.5*pi*(lon(ii)-360)/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(LATraster+33))/(lat(31)-lat(33))))^2);
                end
            end
        elseif mod(LATraster,2) == 0 && LONraster == 32.5
            for ii=LONraster*4-3-2:LONraster*4-3+1
                for jj=LATraster+33-2:LATraster+33+2
                    SOURCEcos(ii,jj)=((cos(0.5*pi*(lon(ii)-lon(LONraster*4-3))/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(LATraster+33))/(lat(31)-lat(33))))^2);
                end
            end
        else
            for ii=LONraster*4-3-2:LONraster*4-3+2
                for jj=LATraster+33-2:LATraster+33+2
                    SOURCEcos(ii,jj)=((cos(0.5*pi*(lon(ii)-lon(LONraster*4-3))/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(LATraster+33))/(lat(31)-lat(33))))^2);
                end
            end
        end
        SOURCEcos=COScoeff*SOURCEcos;
        FORCING = FORCING + SOURCEcos;
    end
   
    if COSforcingsingle==1
        if LONraster == 1
            for ii=LONraster*2-1:LONraster*2-1+2
                for jj=35-2:35+2
                    SOURCEcos(ii,jj)=((cos(0.5*pi*(lon(ii)-lon(LONraster*2-1))/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(35))/(lat(31)-lat(33))))^2);
                end
            end
            for ii=127:128
                for jj=35-2:35+2
                    SOURCEcos(ii,jj)=SOURCEcos(ii,jj)+((cos(0.5*pi*(lon(ii)-360)/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(35))/(lat(31)-lat(33))))^2);
                end
            end
        elseif LONraster == 64
            for ii=LONraster*2-1-2:LONraster*2-1+1
                for jj=35-2:35+2
                    SOURCEcos(ii,jj)=((cos(0.5*pi*(lon(ii)-lon(LONraster*2-1))/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(35))/(lat(31)-lat(33))))^2);
                end
            end
        else
            for ii=LONraster*2-1-2:LONraster*2-1+2
                for jj=35-2:35+2
                    SOURCEcos(ii,jj)=((cos(0.5*pi*(lon(ii)-lon(LONraster*2-1))/(lon(3)-lon(1))))^2)*((cos(0.5*pi*(lat(jj)-lat(35))/(lat(31)-lat(33))))^2);
                end
            end
        end
        SOURCEcos=COScoeff*SOURCEcos;
        FORCING = FORCING + SOURCEcos;
    end
end

if NONDIVadjust==1
[FORCING] = GET_NONDIVforcing(FORCING,VORT_gp_INIT,Fcor,lat);
end

ok=22;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the Laplace equation for vorticity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VORT_t = del2(PSI_t); % Vorticity at (t-1) and (t)

DELTA=set(DELTA,'gp',Z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U_t,V_t] = UVinv(VORT_t,DELTA);
%[U_tm1,V_tm1] = UVinv(VORT_tm1,DELTA);

%======= SAVE OUT INITIAL FIELDS =============

U_gp_t = get(U_t,'gp'); % Grid velocities at (t-1) and (t)
V_gp_t = get(V_t,'gp');

U250=cat(3,U250,U_gp_t);
V250=cat(3,V250,V_gp_t);

PSI_t = del2inv(VORT_t);
VORT_tm1 = del2(PSI_t);

psi250=cat(3,psi250,get(PSI_t,'gp'));

vort250=cat(3,vort250,get(VORT_t,'gp'));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin time stepping  -- Run method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running...')
for nstep = 1:nsteps

U_gp_t = get(U_t,'gp'); % Grid velocities at (t-1) and (t)
V_gp_t = get(V_t,'gp');
    
%=========== ADD CORIOLIS FORCE TO GP VORT AT (t) =========    
    
VORT_gp_t = get(VORT_t,'gp');
ABSVORT_gp_t = VORT_gp_t + Fcor;   

%=========== GET ADVECTION TERMS AT (t) ==============

Uadvect_gp_t = U_gp_t.*ABSVORT_gp_t;
Vadvect_gp_t = V_gp_t.*ABSVORT_gp_t;

Uadvect_t = set(Uadvect_t,'gp',Uadvect_gp_t);
Uadvect_t = shtrana(Uadvect_t);

Vadvect_t = set(Vadvect_t,'gp',Vadvect_gp_t);
Vadvect_t = shtrana(Vadvect_t);

ADVECT_t = div(Uadvect_t,Vadvect_t);

%=========== DO IMPLICIT DAMPING AND TIMESTEP (SPECTRAL) ======

ADVECT_gp_t = get(ADVECT_t,'gp');
VORT_gp_tm1 = get(VORT_tm1,'gp');

if DIFForder==1
    DIFFUSE_t = del2(VORT_t);
elseif DIFForder==2
    DIFFUSE_t = del2(del2(VORT_t));
end
DIFFUSE_gp_t = get(DIFFUSE_t,'gp');

if GETinitTEND==1 && nstep==1
    INITforcing = - (dt)*(ADVECT_gp_t)- dt*DRAGcoeff*VORT_gp_t + dt*DIFFcoeff*DIFFUSE_gp_t + dt*FORCING;
    save(FORCEfname,'INITforcing');
end

if RUNstationary==1
    VORT_gp_tp1 = VORT_gp_t - (dt)*(ADVECT_gp_t)- dt*DRAGcoeff*VORT_gp_t + dt*DIFFcoeff*DIFFUSE_gp_t + dt*FORCING - INITforcing;
   % VORT_tp1 = set(VORT_tp1,'gp',VORT_gp_tp1);
    %VORT_tp1 = shtrana(VORT_tp1);
else 
    VORT_gp_tp1 = VORT_gp_t - (dt)*(ADVECT_gp_t)- dt*DRAGcoeff*VORT_gp_t + dt*DIFFcoeff*DIFFUSE_gp_t;
    %VORT_tp1 = set(VORT_tp1,'gp',VORT_gp_tp1);
    %VORT_tp1 = shtrana(VORT_tp1);
end


%========== DO ROBERT FILTER ON VORT (t) =====================

% VORT_gp_tp1 = get(VORT_tp1,'gp');
% VORT_gp_tp1 = VORT_gp_tp1 - 2*dt*DRAGcoeff*VORT_gp_tp1;
% 
% VORT_gp_t = (1-2*ROBERTcoeff)*VORT_gp_t + ROBERTcoeff*(VORT_gp_tp1 + VORT_gp_tm1);
% 
% VORT_t = set(VORT_t,'gp',VORT_gp_t);
% VORT_t = shtrana(VORT_t);


%========= OUTPUT TO ARRAY/PLOT ==================

%VORT_t = VORT_tp1;
VORT_gp_t = VORT_gp_tp1;
VORT_t = set(VORT_t,'gp',VORT_gp_t);
[U_t,V_t] = UVinv(VORT_t,DELTA);
U_gp_t = get(U_t,'gp'); % Grid velocities at (t-1) and (t)
V_gp_t = get(V_t,'gp');

if (mod(nstep,nout) == 0 ) 
  disp(nstep*dt/86400)
  
  if PLOTfig==1
  
  %figure(1);
  %pcolor(lon,lat,get(VORT_t,'gp')'); shading flat; colorbar; %caxis(1e-5*[-1 1]); %pause(0.2);
  %pause(0.1);
    
  %figure(2);
  %pcolor(lon,lat,get(VORT_t,'gp')' - (squeeze(vort250(:,:,1))')); shading flat; colorbar; %caxis(1e-5*[-1 1]); %pause(0.2);
  %pause(0.1);
  
  %figure(3);
  %pcolor(lon,lat,get(U_t,'gp')'); shading flat; colorbar; %caxis(1e-5*[-1 1]); %pause(0.2);
  %pause(0.1);
  
  %figure(4);
  %pcolor(lon,lat,get(U_t,'gp')' - (squeeze(U250(:,:,1))')); shading flat; colorbar; %caxis(1e-5*[-1 1]); %pause(0.2);
  %pause(0.1);
  
  %figure(5);
  %pcolor(lon,lat,get(V_t,'gp')'); shading flat; colorbar; %caxis(1e-5*[-1 1]); %pause(0.2);
  %pause(0.1);
  
  figure(6);
  pcolor(lon,lat,get(V_t,'gp')' - (squeeze(V250(:,:,1))')); shading flat; colorbar; %caxis(1e-5*[-1 1]); %pause(0.2);
  pause(0.1);
  
  %figure(7);
  %pcolor(lon,lat,FORCING'); shading flat; colorbar; %caxis([-3 3]); %pause(0.2);
  %pause(0.1);
  
  end
    
  U250=cat(3,U250,U_gp_t);
  V250=cat(3,V250,V_gp_t);

  PSI_t = del2inv(VORT_t);
  psi250=cat(3,psi250,get(PSI_t,'gp'));

  vort250=cat(3,vort250,get(VORT_t,'gp'));
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Prognostics are at the new time level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VORT_t = VORT_tp1;
% VORT_tm1 = VORT_t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update diagnostic relations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate velocities by inverting the UV relation
 % [U_t,V_t] = UVinv(VORT_t,DELTA);
  

  
end

% save(['outputs/' OUTfname '.mat'],'vort250','psi250','U250','V250','SOURCE','lat','lon');


vort250=cat(1,vort250(65:128,:,:),vort250(1:64,:,:));
psi250=cat(1,psi250(65:128,:,:),psi250(1:64,:,:));
U250=cat(1,U250(65:128,:,:),U250(1:64,:,:));
V250=cat(1,V250(65:128,:,:),V250(1:64,:,:));
SOURCE=cat(1,SOURCE(65:128,:,:),SOURCE(1:64,:,:));

%=================== WRITE OUTPUT ===============


disp(' Create netCDF4 file')

		% Define axis, for example:
		LON = lon;%-180; % Longitude
        %LON = lon;
		LAT = flipud(lat); % Latitude
		T = 1:size(vort250,3); % Time

		% Create a random field to record in the cdf file dimensions are (Z,Y,X):
% 		C = rand(length(Z),length(LAT),length(LON));

%disp(' File standard set-up')

		% Open file and insert dimensions

		%scope = netcdf.create(['outputs/single/' model '/' MONname '/' OUTfname '.nc'],'netcdf4');
        scope = netcdf.create(['outputs/' OUTfname '.nc'],'netcdf4');
		% Define useful constants:
		NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
		fillValue = -999;

		% Define dimensions:
		dimidY = netcdf.defDim(scope,'latitude',length(LAT));
		dimidX = netcdf.defDim(scope,'longitude',length(LON));

		% Define axis:

		varid = netcdf.defVar(scope,'latitude','double',[dimidY]);
		netcdf.putAtt(scope,varid,'standard_name','latitude');
		netcdf.putAtt(scope,varid,'long_name','Grid latitude');
		netcdf.putAtt(scope,varid,'units','degrees_north');
		netcdf.defVarFill(scope,varid,false,fillValue);
		netcdf.putVar(scope,varid,LAT);

		varid = netcdf.defVar(scope,'longitude','double',[dimidX]);
		netcdf.putAtt(scope,varid,'standard_name','longitude');
		netcdf.putAtt(scope,varid,'long_name','Grid longitude');
		netcdf.putAtt(scope,varid,'units','degrees_east');
		netcdf.defVarFill(scope,varid,false,fillValue);
		netcdf.putVar(scope,varid,LON);


        	dimidT = netcdf.defDim(scope,'time',length(T));

        varid = netcdf.defVar(scope,'time','double',[dimidT]);			
        netcdf.putAtt(scope,varid,'standard_name','time');
        netcdf.putAtt(scope,varid,'long_name','Time of measurements');
        netcdf.putAtt(scope,varid,'units','days since start');
        netcdf.defVarFill(scope,varid,false,fillValue);
        netcdf.putVar(scope,varid,T);

        
	% Insert global attributes

		netcdf.putAtt(scope,NC_GLOBAL,'title','BT vorticity equation output')
		netcdf.putAtt(scope,NC_GLOBAL,'long_title','')
		netcdf.putAtt(scope,NC_GLOBAL,'comments','Should be ok...')
		netcdf.putAtt(scope,NC_GLOBAL,'institution','Oxford')
		netcdf.putAtt(scope,NC_GLOBAL,'source','My barotropic model')

		netcdf.putAtt(scope,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'))
		netcdf.putAtt(scope,NC_GLOBAL,'CreatedBy',getenv('LOGNAME'))
		netcdf.putAtt(scope,NC_GLOBAL,'MatlabSource',which('Matlab_netcdf_update_sample'))
        
        

        varname = 'V';
        long_name = 'Meridional velocity (250hPa)';
        unit = 'm/s';

        varid = netcdf.defVar(scope,varname,'double',[dimidX,dimidY,dimidT]);			
        netcdf.putAtt(scope,varid,'long_name',long_name);
        netcdf.putAtt(scope,varid,'units',unit);
        netcdf.defVarFill(scope,varid,false,fillValue);
        v = flip(circshift(V250,-64,1),2);
		v(isnan(v)) = fillValue;
		netcdf.putVar(scope,varid,v);
        
        varname = 'U';
        long_name = 'Zonal velocity (250hPa)';
        unit = 'm/s';

        varid = netcdf.defVar(scope,varname,'double',[dimidX,dimidY,dimidT]);			
        netcdf.putAtt(scope,varid,'long_name',long_name);
        netcdf.putAtt(scope,varid,'units',unit);
        netcdf.defVarFill(scope,varid,false,fillValue);
        v = flip(circshift(U250,-64,1),2);
		v(isnan(v)) = fillValue;
		netcdf.putVar(scope,varid,v);
        
        varname = 'PSI';
        long_name = 'Streamfunction (250hPa)';
        unit = 'm^2/s';

        varid = netcdf.defVar(scope,varname,'double',[dimidX,dimidY,dimidT]);			
        netcdf.putAtt(scope,varid,'long_name',long_name);
        netcdf.putAtt(scope,varid,'units',unit);
        netcdf.defVarFill(scope,varid,false,fillValue);
        v = flip(circshift(psi250,-64,1),2);
		v(isnan(v)) = fillValue;
		netcdf.putVar(scope,varid,v);
        
        varname = 'VORTICITY';
        long_name = 'Vorticity (250hPa)';
        unit = '1/s';

        varid = netcdf.defVar(scope,varname,'double',[dimidX,dimidY,dimidT]);			
        netcdf.putAtt(scope,varid,'long_name',long_name);
        netcdf.putAtt(scope,varid,'units',unit);
        netcdf.defVarFill(scope,varid,false,fillValue);
        v = flip(circshift(vort250,-64,1),2);
		v(isnan(v)) = fillValue;
		netcdf.putVar(scope,varid,v);
        
        varname = 'FORCING';
        long_name = 'Rossby wave source forcing (250hPa)';
        unit = '1/s^2';

        varid = netcdf.defVar(scope,varname,'double',[dimidX,dimidY]);			
        netcdf.putAtt(scope,varid,'long_name',long_name);
        netcdf.putAtt(scope,varid,'units',unit);
        netcdf.defVarFill(scope,varid,false,fillValue);
        v = flip(circshift(SOURCE,-64,1),2);
		v(isnan(v)) = fillValue;
		netcdf.putVar(scope,varid,v);
        
        
	% Check what we've just donewithou
		netcdf.close(scope)
        
        %unix(['cp ../Vars/' OUTfname '.nc ../../../../MCAcode/Vars/']);


ok=22;

end


