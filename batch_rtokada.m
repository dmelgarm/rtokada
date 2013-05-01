function batch_rtokada

%Diego Melgar 01/2013

%set up paralell run
%matlabpool
%Set barch run
stations='stations_kal.xy';
%runName='RTOkada_HiResW__checker_kal_LCopt'
runName='RTOkada_HiResW_kal_gps'
%stasuffix='raw'
%stasuffix='checker_kal';
stasuffix='kal_gps';
%stasuffix='checker'
workpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada'; 
outpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada/output/';
%lambda=logspace(-0.30103,3.6989,200);
%lambda=logspace(-1,2,200);
%lambda=logspace(0,log10(500),200);
lambda=4.7;
weightflag=1;
N=length(lambda);
cd(workpath)
% load('green.mat')
% load('greenSF.mat')
% load('green_small.mat')
load('green_small_kal.mat')
load('greenSF_small.mat')
Gp=G;
Gs=GSF;

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
    [l(k) L2(k) LS(k) Mo(k) Mw(k) VR(k) GCV(k) ABIC(k)]=rtokada(workpath,outpath,runName,runID{k},stasuffix,lambda(k),Gp,Gs,weightflag,stations);
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