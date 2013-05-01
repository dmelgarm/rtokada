function wetgreen

%DMM 04/2013
%This code will make the seafloor displacements for 1m of strike slip and
%1m of dip slip motion on all subfaults. These dispalcements will then be
%used in GeoClaw to compute kernels at all gauges.

workpath='~/Research/Data/Tohoku/RTOkada'
outpath=[workpath '/wetgreen/'];
runName='tohoku_2patch_0010_0037';
%Make slip
slip=1; %in m
%Load Seafloor GFs
cd(workpath)
load greenSF_small.mat
[latsf lonsf]=textread('seafloor.xy','%f%f');
%Make displacements
npatch=size(GSF,2);
fault=zeros(npatch,1);
%Do strike slip first
for np=1:2:npatch
    NP=floor(np/2)+1
    if NP<10
        subfault=['000' num2str(NP)];
    elseif NP<100
        subfault=['00' num2str(NP)];
    elseif NP<1000
        subfault=['0' num2str(NP)];
    else
        subfault=num2str(NP);
    end
    fid=fopen([outpath runName '.' subfault '.ss.wgf'],'wt');
    fault(np)=fault(np)+slip;
    D=GSF*fault;
    fault(np)=fault(np)-slip;
    %Write  GeoClaw dtopo type 3 file
    for i=1:length(lonsf)
        fprintf(fid,'%1.0f %1.6f %1.6f %1.12f\n',0,lonsf(i),latsf(i),0);
    end
    nnn=1;
    for i=1:length(lonsf)
        fprintf(fid,'%1.0f %1.6f %1.6f %1.12f\n',1,lonsf(i),latsf(i),D((nnn-1)*3+3));
        nnn = nnn+1;
    end
    fclose(fid);
end
%Now dip-slip
for np=2:2:npatch
%for np=20;
    NP=floor(np/2);
    if NP<10
        subfault=['000' num2str(NP)];
    elseif NP<100
        subfault=['00' num2str(NP)];
    elseif NP<1000
        subfault=['0' num2str(NP)];
    else
        subfault=num2str(NP);
    end
    fid=fopen([outpath runName '.' subfault '.ds.wgf'],'wt');
    fault(np)=fault(np)+slip;
    %fault(np+54)=fault(np+54)+slip;
    D=GSF*fault;
    fault(np)=fault(np)-slip;
    %Write  GeoClaw dtopo type 3 file
    for i=1:length(lonsf)
        fprintf(fid,'%1.0f %1.6f %1.6f %1.12f\n',0,lonsf(i),latsf(i),0);
    end
    nnn=1;
    for i=1:length(lonsf)
        fprintf(fid,'%1.0f %1.6f %1.6f %1.12f\n',1,lonsf(i),latsf(i),D((nnn-1)*3+3));
        nnn = nnn+1;
    end
    fclose(fid);
end

