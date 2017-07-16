clear
clear all

writetofile=true;
plotfigure=true;

% Resolution
% nlev stands for number of vertical levels
% nlat stands for horizontal resolution, 32 = T21, 64 = T42, 128 = T85, 256 = T170 
nlev=30;
nlat=32;                %note: nlat<1000

% off equator heat parameter ?unit degree?
phi00=0;  %233phi0tag

ps_const=1011.00;       % surface pressure unit hpa
p0=1000.0;              % reference level

% Do not change following code unless you really know what you are doing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up grid point in gaussian gird 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi0=pi/180*phi00; %convert to rad unit
nlon=nlat*2;
nugp=nlat*nlon;
sigma=0.5/nlev:1/nlev:1-0.5/nlev;   % sigma level 
plev=sigma*ps_const;                % pressure level 

[xlat,dlat,sinc]=gauss2lats(nlat); % gaussian grid
lat=[xlat -fliplr(xlat)]; % -90:90 gaussian latitude
phi=lat/180*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. set up radiative restoration temperature 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ta1=zeros(nlon,nlat,nlev);  % constant 
Ta2=zeros(nlon,nlat,nlev);  % time variable for seasonal cycle

% Walker and Schneider 2005 (identical to Held Suarze 1994, except that
% the off equator heating parameter phi0
for k=1:nlev
    for j=1:nlat
        Ta1(:,j,k)=max(200,(315-60*(sin(phi(j))^2-2*sin(phi0)*sin(phi(j)))-10*log(plev(k)/p0)*cos(phi(j)-phi0)^2)*(plev(k)/p0)^(2/7));
    end
end

% % PKU_GCM aqua set for temperature (no idea...)
% Tre=?;

% % Schneider_Bordoni_2008 & Bordoni_Schneider_2010
% Tres=max(200,350-112.5*(sin(phi).^2-2*sin(phi0)*sin(phi)));

% % Walker_Schneider_2006 
% Tres=260+120*(cos(phi).^2+2*sin(phi0)*(1+sin(phi)));
% for k=1:nlev
%     for j=1:nlat
%         Ta1(:,j,k)=200*(1+((Tres(j)/200)^4-1)*(plev(k)/p0)^3.5)^(1/4);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. set up temperature restoration time scale -- taur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taurlll=zeros(nlon,nlat,nlev);
for k=1:nlev

    % Held Suarze 1994
    for j=1:nlat
        taurlll(:,j,k)=86400/(1/40+(1/4-1/40)*max(0,(sigma(k)-0.7)/(1-0.7))*cos(phi(j))^4);
    end
    
%     %set all to 50 days like paper of Walker and Schneider 2006 
%     taurlll(:,:,k)=50*86400;
%     %Original set in PKU_GCM(puma)
%     taurlll(:,:,k)=min(30,50*atan(1-sigma(k)))*86400;

%     %Schneider(2004)
%     for j=1:nlat
%         taurlll(:,j,k)=86400/(1/50+(1/7-1/50)*max(0,(sigma(k)-0.85)/(1-0.85))*cos(phi(j))^8);
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. set up topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
topo=zeros(nlon,nlat,nlev);  % set topography for auqa experiment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. set up initial surface pressure in order to reduce model relax time due to topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps=ones(nlon,nlat,nlev);
ps=ps*ps_const;   % 101100.0 is default value of PSURF_EARTH in PKU_GCM

%% plot section
if (plotfigure == true)
    [X,Y]=meshgrid(lat,flip(sigma));
    figure(1)
    subplot(3,2,1)
    plot(squeeze(Ta1(1,nlat/2,:)),sigma)
    title('Treq vertical profile in equator')
    set(gca, 'YDir','reverse');
    ylabel('sigma')
 %%  
    sb1=subplot(3,2,2);
    contourf(X,Y,flipud(squeeze(Ta1(1,:,:))'),12)
    set(gca, 'YDir','reverse');
    colormap(sb1,jet(12))
    title('Constant radiative temperature(K)')
    caxis([200 320])
    colorbar;
%%    
    pt=squeeze(Ta1(1,:,:));
    for k=1:nlev
        pt(:,k)=pt(:,k).*(plev(k)/p0)^(-2/7);
    end
    sb2=subplot(3,2,[3 4]);
    contourf(X,Y,flipud(pt(:,:)'),24)
    title('Constant radiative potential temperature(K)')
    set(gca, 'YDir','reverse');
    %caxis([315 320])
    colormap(sb2,jet(12))
    colorbar;
    
    subplot(3,2,5)
    plot(squeeze(taurlll(1,nlat/2,:)/86400),sigma)
    title('Taur vertical profile in equator')
    set(gca, 'YDir','reverse');
    ylabel('sigma')
    
    sb3=subplot(3,2,6);
    contourf(X,Y,flipud(squeeze(taurlll(1,:,:)/86400)'),12)
    colormap(sb3,jet(12));
    title('Relaxaton time(day)')
    set(gca, 'YDir','reverse');
    colorbar
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write all to file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (writetofile == true) 
    texttmp='Write all to files...'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write into file Nxxx_surf_0121.sra
    % constant radiative restoration temperature -- gr1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kcode=121;                          % kcode=121 stand for restoration temperature
    filename=[  'N'...
                num2str(sprintf('%0.3i',nlat))...
                '_surf_'...
                num2str(sprintf('%0.4i',kcode))...
                '.sra'...
                num2str(sprintf('%0.2i',phi00))]
    fileID=fopen(filename,'w');               
    %head in file
    ihead=zeros(1,8);
    ihead(1)=kcode;                     % kcode=121 stand for restoration temperature
    ihead(3)=str2num(datestr(now,'yyyymmdd'));   % YearMonthDay
    ihead(4)=str2num(datestr(now,'HHMM'));       % Hour Minites
    ihead(5)=nlon;                      % 1.dimension
    ihead(6)=nlat;                      % 2.dimension
    ihead(7)=1;
    ihead(8)=0;
    for k=1:nlev
        ihead(2)=k;                     %level
        fprintf(fileID,'%10i%10i%10i%10i%10i%10i%10i%10i\n',ihead);
        fprintf(fileID,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',reshape(Ta1(:,:,k),[],1));
    end
    fclose(fileID);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write into file Nxxx_surf_0122.sra 
    % variable radiative restoration temperature -- gr2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kcode=122;                          % kcode=122 stand for seasonal part of restoration temperature
    filename=[  'N'...
                num2str(sprintf('%0.3i',nlat))...
                '_surf_'...
                num2str(sprintf('%0.4i',kcode))...
                '.sra']
    fileID=fopen(filename,'w');              
    %head in file
    ihead=zeros(1,8);
    ihead(1)=kcode;                     
    ihead(3)=str2num(datestr(now,'yyyymmdd'));   % YearMonthDay
    ihead(4)=str2num(datestr(now,'HHMM'));       % Hour Minites
    ihead(5)=nlon;                      % 1.dimension
    ihead(6)=nlat;                      % 2.dimension
    ihead(7)=1;
    ihead(8)=0;
    for k=1:nlev
        ihead(2)=k;                     %level
        fprintf(fileID,'%10i%10i%10i%10i%10i%10i%10i%10i\n',ihead);
        fprintf(fileID,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',reshape(Ta2(:,:,k),[],1));
    end
    fclose(fileID);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write into file Nxxx_surf_0123.sra 
    % temperature restoration time scale -- taurlll
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kcode=123;                          % kcode=123 stand for temperature restoration time scale
    filename=[  'N'...
                num2str(sprintf('%0.3i',nlat))...
                '_surf_'...
                num2str(sprintf('%0.4i',kcode))...
                '.sra']
    fileID=fopen(filename,'w');     
    %head in file
    ihead=zeros(1,8);
    ihead(1)=kcode;                     
    ihead(3)=str2num(datestr(now,'yyyymmdd'));   % YearMonthDay
    ihead(4)=str2num(datestr(now,'HHMM'));       % Hour Minites
    ihead(5)=nlon;                      % 1.dimension
    ihead(6)=nlat;                      % 2.dimension
    ihead(7)=1;
    ihead(8)=0;
    for k=1:nlev
        ihead(2)=k;                     %level
        fprintf(fileID,'%10i%10i%10i%10i%10i%10i%10i%10i\n',ihead);
        fprintf(fileID,'%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n',reshape(taurlll(:,:,k),[],1));
    end
    fclose(fileID);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write into file Nxxx_surf_0129.sra 
    % Topography -- topo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kcode=129;                          % kcode=129 stand for Topography
    filename=[  'N'...
                num2str(sprintf('%0.3i',nlat))...
                '_surf_'...
                num2str(sprintf('%0.4i',kcode))...
                '.sra']
    fileID=fopen(filename,'w');     
    %head in file
    ihead=zeros(1,8);
    ihead(1)=kcode;                     
    ihead(3)=str2num(datestr(now,'yyyymmdd'));   % YearMonthDay
    ihead(4)=str2num(datestr(now,'HHMM'));       % Hour Minites
    ihead(5)=nlon;                      % 1.dimension
    ihead(6)=nlat;                      % 2.dimension
    ihead(7)=1;
    ihead(8)=0;
    for k=1:nlev
        ihead(2)=k;                     %level
        fprintf(fileID,'%10i%10i%10i%10i%10i%10i%10i%10i\n',ihead);
        fprintf(fileID,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',reshape(topo(:,:,k),[],1));
    end
    fclose(fileID);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write into file Nxxx_surf_0134.sra 
    % initial surface pressure -- taur
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kcode=134;                          % kcode=134 stand for initial surface pressure
    filename=[  'N'...
                num2str(sprintf('%0.3i',nlat))...
                '_surf_'...
                num2str(sprintf('%0.4i',kcode))...
                '.sra']
    fileID=fopen(filename,'w');     
    %head in file
    ihead=zeros(1,8);
    ihead(1)=kcode;                     
    ihead(3)=str2num(datestr(now,'yyyymmdd'));   % YearMonthDay
    ihead(4)=str2num(datestr(now,'HHMM'));       % Hour Minites
    ihead(5)=nlon;                      % 1.dimension
    ihead(6)=nlat;                      % 2.dimension
    ihead(7)=1;
    ihead(8)=0;
    for k=1:nlev
        ihead(2)=k;                     %level
        fprintf(fileID,'%10i%10i%10i%10i%10i%10i%10i%10i\n',ihead);
        fprintf(fileID,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',reshape(ps(:,:,k),[],1));
    end
    fclose(fileID);
end

%%
% This function provides latitudes on a Gaussian grid from the
% number of latitude lines.
%
% Source: https://www.mathworks.com/matlabcentral/fileexchange/2586-gauss2lats
% the calling format is:
%               [xlat,dlat,sinc]=gauss2lats(nlat)
%     where: input is nlat = the number of latitude lines, e.g 94 in a T62 grid
%                             with 192 longitudes and 94 latitudes
%            outputs are:
%                      xlat = the latitudes
%                      dlat = latitude spacing
%                      sinc = sine of colatitudes (cosine of latitudes)
% 
% Adapted for Matlab from the NCAR Fortran program by Tom Holt, 23/10/2002.
%
% Possible problems:
%          1) The output may appear slightly different from other
%          estimates. I am indebted to Lee Panetta who provided the
%          following test results (21.10.2003), comparing Gauss2lats on a
%          G4 Mac with another algorithm on a Cray J90:
%                   agreement to 15 decimal places except at endpoints
%                   where agreement was to 13 decimal places.
%          2) The routine may not converge to a solution on some grids.

function [xlat,dlat,sinc]=gauss2lats(nlat)
acon=180.0/pi;
% convergence criterion for iteration of cos latitude
xlim=1.0e-7;
% initialise arrays
for i=1:720
    cosc(i)=0.;
    gwt(i)=0.;
    sinc(i)=0.;
    colat(i)=0.;
    wos2(i)=0.;
end
% the number of zeros between pole and equator
nzero=nlat/2;
% set first guess for cos(colat)
for i=1:nzero;
    cosc(i)=sin((i-0.5)*pi/nlat+pi*0.5);
end
% constants for determining the derivative of the polynomial
fi=nlat;
fi1=fi+1.0;
a=fi*fi1/sqrt(4.0*fi1*fi1-1.0);
b=fi1*fi/sqrt(4.0*fi*fi-1.0);
%loop over latitudes, iterating the search for each root
for i=1:nzero
    % determine the value of the ordinary Legendre polynomial for the current guess root
    g=gord(nlat,cosc(i));
    % determine the derivative of the polynomial at this point
    gm=gord(nlat-1,cosc(i));
    gp=gord(nlat+1,cosc(i));
    gt=(cosc(i)*cosc(i)-1.0)/(a*gp-b*gm);
    % update the estimate of the root
    delta=g*gt;
    cosc(i)=cosc(i)-delta;
 
    % if convergence criterion has not been met, keep trying
    while abs(delta) > xlim
        g=gord(nlat,cosc(i));
        gm=gord(nlat-1,cosc(i));
        gp=gord(nlat+1,cosc(i));
        gt=(cosc(i)*cosc(i)-1.0)/(a*gp-b*gm);
        delta=g*gt;
        cosc(i)=cosc(i)-delta;   
    end
    % determine the Gaussian weights
    c=2.0*(1.0-cosc(i)*cosc(i));
    d=gord(nlat-1,cosc(i));
    d=d*d*fi*fi;
    gwt(i)=c*(fi-0.5)/d;
end
% determine the colatitudes and sin(colat) and weights over sin**2
for i=1:nzero
    colat(i)=acos(cosc(i));
    sinc(i)=sin(colat(i));
    wos2(i)=gwt(i)/(sinc(i)*sinc(i));
end
% if nlat is odd, set values at the equator
if mod(nlat,2) ~= 0 
    i=nzero+1;
    cosc(i)=0.0;
    c=2.0;
    d=gord(nlat-1,cosc(i));
    d=d*d*fi*fi;
    gwt(i)=c*(fi-0.5)/d;
    colat(i)=pi*0.5;
    sinc(i)=1.0;
    wos2(i)=gwt(i);
end
% determine the southern hemisphere values by symmetry
for i=nlat-nzero+1:nlat
    cosc(i)=-cosc(nlat+1-i);
    gwt(i)=gwt(nlat+1-i);
    colat(i)=pi-colat(nlat+1-i);
    sinc(i)=sinc(nlat+1-i);
    wos2(i)=wos2(nlat+1-i);
end
ylat=90.;
% calculate latitudes and latitude spacing
for i=1:nzero
    xlat(i)=acos(sinc(i))*acon;
    dlat(i)=xlat(i)-ylat;
    ylat=xlat(i);
end
end

% This function calculates the value of an ordinary Legendre polynomial at a latitude.
%          inputs are:
%                n = the degree of the polynomial
%                x = cos(colatitude)
%         outputs are:
%              ggg = the value of the Legendre polynomial of degree n at 
%                    latitude asin(x)
%
function [ggg]=gord(n,x)
% determine the colatitude
colat=acos(x);

c1=sqrt(2.0);

for i=1:n
    c1=c1*sqrt(1.0-1.0/(4*i*i));
end
fn=n;
ang=fn*colat;
s1=0.0;
c4=1.0;
a=-1.0;
b=0.0;
for k=0:2:n
    if k==n 
        c4=0.5*c4;
    end
    s1=s1+c4*cos(ang);
    a=a+2.0;
    b=b+1.0;
    fk=k;
    ang=colat*(fn-fk-2.0);
    c4=(a*(fn-b+1.0)/(b*(fn+fn-a)))*c4;
end
ggg=s1*c1;
end