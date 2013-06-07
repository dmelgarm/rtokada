function batch_rtokada

%Diego Melgar 01/2013

%set up paralell run
%matlabpool
%Set barch run
stations='stations_kal.xy';
runName='RTOkada_wv';
%runName='RTOkada_HiResW__checker_kal_LCopt'
%runName='RTOkada_HiResW_kal_weightgps'
%stasuffix='raw'
%stasuffix='checker_kal';
%stasuffix='kal_weightgps';
stasuffix='kal';
%stasuffix='checker'
workpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada'; 
outpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada/output/';
%lambda=logspace(-0.30103,3.6989,200);
lambda=logspace(-1.30103,0.69897,200);
%lambda=logspace(-1,2,200);
%lambda=logspace(0,log10(500),200);
waveflag=1;   %Use wave gauges
weightflag=1;
N=length(lambda);
cd(workpath)
% load('green.mat')
% load('greenSF.mat')
% load('green_small.mat')
load('green_small_kal.mat')
load('greenSF_small.mat')
load('tohoku_wvGF_60min.mat')
Gp=G;
Gs=GSF;
Gw=Gwv;

numsta=size(Gp,2)/3;
%Get fault to output point distances
lono = 143.05;%starting longitude (corresponding to x=0) - if using large latitude extent, set these to center of area
lato = 37.5;%starting latitude (corresponding to y=0)
[llat,llon] = degreelen(lato);%lengths of degree of lat and lon
[site latinv,loninv]=textread(stations,'%f %f %f');
%Fault to station
% [xs,ys,zs,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=textread('small_fault.dat','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
% for i = 1:length(xs)
%     for j = 1:numsta
%         xrs(i,j) = (loninv(j)-xs(i))*llon;%x distance, in m
%         yrs(i,j) = (latinv(j)-ys(i))*llat;%y distance, in m
%         zrs(i,j) = -zs(i);%z distance, in m
%     end
% end
% % Fault to seafloor ditances
% [latsf lonsf]=textread('seafloor.xy','%f%f');
% for i = 1:length(xs)
%     for j = 1:length(latsf)
%         xrs(i,j) = (lonsf(j)-xs(i))*llon;%x distance, in m
%         yrs(i,j) = (latsf(j)-ys(i))*llat;%y distance, in m
%         zrs(i,j) = -zs(i);%z distance, in m
%     end
% end

%batch run
for k=1:N
    if k<10
        runID{k}=['000' num2str(k)];
    elseif k<100
        runID{k}=['00' num2str(k)];
    elseif k<1000
        runID{k}=['0' num2str(k)];
    else
        runID=num2str(k);
    end
end
tic
for k=1:N
    k
    %tic
    if waveflag==0 %Only dispalcement ivnersion
        [l(k) L2(k) LS(k) Mo(k) Mw(k) VR(k) GCV(k) ABIC(k)]=rtokada(workpath,outpath,runName,runID{k},stasuffix,lambda(k),Gp,Gs,weightflag,stations);
    else %Also ivnert wave gauges
        [l(k) L2(k) LS(k) Mo(k) Mw(k) VR(k) GCV(k) ABIC(k)]=rtokadawv(workpath,outpath,runName,runID{k},stasuffix,lambda(k),Gp,Gs,Gw,weightflag,stations,tGF,gauges);
    end
    %toc
end
toc
for k=1:N
    %Write run log
    fid=fopen([outpath runName '.' runID{k} '.log'],'wt');
    fprintf(fid,'%s %s\n','runID',runID{k});
    fprintf(fid,'%s %1.5f\n','Smoothing', lambda(k));
    fprintf(fid,'%s %1.5f\n','Misfit-L2', L2(k));
    fprintf(fid,'%s %1.5f\n','Solution-Semi-Norm', LS(k));
    fprintf(fid,'%s %1.5f\n','Moment(Nm)', Mo(k));
    fprintf(fid,'%s %1.5f\n','Moment-Magnitude', Mw(k));
    fprintf(fid,'%s %1.5f\n','VR(%)', VR(k));
    fprintf(fid,'%s %1.5f\n','GCV', GCV(k));
    fprintf(fid,'%s %1.5f\n','ABIC', ABIC(k));
    fclose(fid);
end
%matlabpool close