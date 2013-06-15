function batch_rtokada(min)

%Diego Melgar 01/2013
%
% Driver script to run rtokadawv() multiple times and sort out output. To run you must define the parameters 
% contained within the lines defined by "~~~" they are as follows:
%
% RUN PARAMETERS
%
% stations - String containing the filename of the lat,lon coordinates of the coseismic offset locations
% runName - A string identifier for the run
% stasuffix - A string with a suffix modifier for the individuals tation file names. This is useful if you have a station 
%             with two instances (each with it's own file), where the offsets are the same but the weights on each 
%             component change.
% lambda - Vector containing the smoothing parameters to be used for the batch run.
% waveflag - 1 to use wav gauge measurements, 0 to ignore.
% coseisflag - 1 to use coseismic offsets, 0 to ignore
% Gcoseis - Load coseismic offset Green functions
% Gseafl -  Load seafloor uplift Green functions
% Gwave - Load wave gauge Green functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
workpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada';   %Where do the files sued in the inversion live
cd(workpath)
outpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada/output/';  %Output directory 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~  JOINT FULL GPS - WAVE GAUGE RUN  ~~~~~~~~~~~~~~~~~~~~~~~
% stations='stations.xy';
% runName='RTOkada_gpswv';
% stasuffix='raw';
% lambda=logspace(-2,1,200);
% lambda=0.22;
% waveflag=1;   %Use wave gauges
% coseisflag=1; %Use coseismic offsets
% Gcoseis='green_small.mat';
% Gseafl='greenSF_small.mat';
% Gwave='tohoku_wvGF_60min.mat';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~ JOINT KALMAN - WAVE GAUGE RUN   ~~~~~~~~~~~~~~~~~~~~~
stations='stations_kal.xy';
runName=['RTOkada_kalwv' min];
stasuffix='kal';
lambda=logspace(log10(0.005),log10(5),200);
waveflag=1;   %Use wave gauges
coseisflag=1; %Use coseismic offsets
Gcoseis='green_small_kal.mat';
Gseafl='greenSF_small.mat';
Gwave=['tohoku_wvGFnoSS_' min '.mat'];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~   FULL GPS NETWORK RUN   ~~~~~~~~~~~~~~~~~~~~~~~~~
% stations='stations.xy'; 
% runName='RTOkada_allgps_opt';
% stasuffix='raw';
% lambda=logspace(log10(0.001),log10(1),200);
% lambda=0.09;
% waveflag=0;   %Use wave gauges
% coseisflag=1; %Use coseismic offsets
% Gcoseis='green_small.mat';
% Gseafl='greenSF_small.mat';
% Gwave='tohoku_wvGF_60min.mat';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~  KALMAN COSEISMICS RUN  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% stations='stations_kal.xy';
% runName='RTOkada_kal_opt';
% stasuffix='kal';
% lambda=logspace(log10(0.0001),0,200);
% lambda=0.01;
% waveflag=0;   %Use wave gauges
% coseisflag=1; %Use coseismic offsets
% Gcoseis='green_small_kal.mat';
% Gseafl='greenSF_small.mat';
% Gwave='tohoku_wvGF_60min.mat';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~  SINGLE RUNS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% stations='stations_kal.xy';
% runName='RTOkada_test';
% stasuffix='kal';
% lambda=logspace(log10(0.0001),0,200);
% lambda=0.219;
% waveflag=1;   %Use wave gauges
% coseisflag=1; %Use coseismic offsets
% Gcoseis='green_small_kal.mat';
% Gseafl='greenSF_small.mat';
% Gwave='tohoku_wvGF_60min.mat';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% GO!
load(Gcoseis)
load(Gseafl)
load(Gwave)
Gp=G;
Gs=GSF;
Gw=Gwv;
N=length(lambda);
cd(workpath)
numsta=size(Gp,2)/3;
%batch run
for k=1:N
    display(['Working on inversion ' num2str(k) ' of ' num2str(N) '...'])
    if k<10
        runID{k}=['000' num2str(k)];
    elseif k<100
        runID{k}=['00' num2str(k)];
    elseif k<1000
        runID{k}=['0' num2str(k)];
    else
        runID=num2str(k);
    end
    dataflag=[coseisflag waveflag];
    [l(k) L2(k) LS(k) Mo(k) Mw(k) GCV(k) AIC(k) VRgps(k) rmswv(k)]=rtokadawv(workpath,outpath,runName,runID{k},stasuffix,lambda(k),Gp,Gs,Gw,stations,tGF,gauges,dataflag);
    %Write run log
    fid=fopen([outpath runName '.' runID{k} '.log'],'wt');
    fprintf(fid,'%s %s\n','runID',runID{k});
    fprintf(fid,'%s %1.5f\n','Smoothing', lambda(k));
    fprintf(fid,'%s %1.5f\n','Misfit-L2', L2(k));
    fprintf(fid,'%s %1.5f\n','Solution-Semi-Norm', LS(k));
    fprintf(fid,'%s %1.5f\n','Moment(Nm)', Mo(k));
    fprintf(fid,'%s %1.5f\n','Moment-Magnitude', Mw(k));
    fprintf(fid,'%s %1.5f\n','VRgps', VRgps(k));
    fprintf(fid,'%s %1.5f\n','RMSwv', rmswv(k));
    fprintf(fid,'%s %1.5f\n','AIC', AIC(k));
    fprintf(fid,'%s %1.5f\n','GCV', GCV(k));
    fclose(fid);
end