function coseisgreen

% D. Melgar. (06/2013)
%
% Function to get Green's functions for coseismic offsets. This is reworked
% from code originally inside B.Crowell's rtokada routine. It computes
% N,E,U for every subfault-station pair. It also computes the GFs for every
% seafloor dispalcement point - subfault pair.
%
% RUN PARAMETERS
% Gcoseis - Name of output file for coseismic GFs
% Gseaf - Name of output file for seafloor GFs
% fault - Name of inpout file containing the subfault metadata
% fault - Name of inpout file containing the subfault coordinates
% stations - Name of file containing the coseismic station coordinates
% seafloor - Name of file containing the lat lons of seafloor displacement points
% lono - Origin longitude (x=0)
% lato - Origin latitude (y=0)


%~~~~~~~~~~~~~~~     RUN PARAMETERS    ~~~~~~~~~~~~~~~~~
Gcoseis='green_small_kal.mat';
Gseafl='greenSF_small_kal.mat';
faultdef='faults_def_small.txt';
fault='small_fault.dat';
stations='stations_kal.xy'
seafloor='seafloor.xy';
lono = 143.05;
lato = 37.5;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



[llat,llon] = degreelen(lato);%lengths of degree of lat and lon
[site latinv,loninv]=textread(stations,'%f %f %f');
f=load(faultdef); %No of fault elements
ast=f(1);%along strike elements
adi=f(2);%along dip elements

%Load fault file, its format is:
%Coordiantes of the centre of the patch, coordinates of the 4 corners,strike,dip,along strike elngth,along dip length, total area
[xs,ys,zs,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=textread(fault,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

%Compute distances
for i = 1:length(xs)
    for j = 1:length(site)
        xrs(i,j) = (loninv(j)-xs(i))*llon;%x distance, in m
        yrs(i,j) = (latinv(j)-ys(i))*llat;%y distance, in m
        zrs(i,j) = -zs(i);%z distance, in m
    end
end
%Compute GFs for coseismic offsets
[G] = okada_green(xrs,yrs,zrs,strike,dip,width.*1000,len.*1000);
save(Gcoseis,'G');

%Now GF for seafloor dispalcements
[latsf lonsf]=textread(seafloor,'%f%f');
for i = 1:length(xs)
    for j = 1:length(latsf)
        xrs(i,j) = (lonsf(j)-xs(i))*llon;%x distance, in m
        yrs(i,j) = (latsf(j)-ys(i))*llat;%y distance, in m
        zrs(i,j) = -zs(i);%z distance, in m
    end
end
[GSF] = okada_green(xrs,yrs,zrs,strike,dip,width.*1000,len.*1000);
save(Gseafl,'GSF');