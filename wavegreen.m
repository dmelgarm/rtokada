function wavegreen

%DMM 05/2013
%
%Assemble matrix of wave gauge green functions

%Run parameters
dt=15; %Sampling rate ins econds
lt=3600; %length of GFs in seconds
Nrecord=(3600/dt)+1; %length of time series
Nf=189;%Number of fault planes
Nsta=15; %Ignoring BY.MIC
dir='/Volumes/Elements/Tsunami/GFs/';  %maind irectory
suffix='GF';
sta_suffix='_60min';
savef='~/Research/Data/Tohoku/RTOkada/tohoku_wvGF_60min';
datapath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada/gauges/';
%Load gauges table (GeoClaw gauge codes vs. actual gauge codes)
[numg a a nameg a a]=textread('/Users/dmelgarm/Research/Data/Tohoku/Gauges/gauges.dat','%f%f%f%s%s%s');
tGF=[];
%Retrieve them
for k=1:Nf
    display(['Assembling subfault ' num2str(k) ' of ' num2str(Nf) '...'])
    if k<10
        subfault=['000' num2str(k)];
    elseif k<100
        subfault=['00' num2str(k)];
    elseif k<1000
        subfault=['0' num2str(k)];
    else
        subfault=num2str(k);
    end
    %Read strike slip
    [sta amr a1 a2 a3 a4 a5]=textread([dir suffix 'ss_' subfault '/_output/fort.gauge'],'%f%f%f%f%f%f%f');
    staG=unique(sta);
    staG=setxor(staG,1016); %Ignore BY.MIC
    Gsta=[];
    for k2=1:Nsta %Loop over all stations for that subfault
        i=find(sta==staG(k2)); %Get next stations data
        t=a1(i);
        eta=(a5(i));
        tG=0:dt:lt;
        etaG=interp1(t,eta,tG);
        %Now equalize to what teh data for that station actually is
        j=find(numg==staG(k2)); %What is the wave gauge code?
        [tdata etadata]=textread([datapath nameg{j} sta_suffix '.txt'],'%f%f'); %load actual data
        tdata=tdata(2:end); %First sample is inversion weight
        etaG=interp1(tG,etaG,tdata); %Reinterpolate to what is actually int he data
        %
        Gsta=[Gsta ; etaG]; %Append current station
        if k==1 %Save time vector for GFs
            tGF=[tGF ; tdata];
        end
    end
    G(:,2*k-1)=Gsta; %Add to main matrix
    
    %Read dip slip
    [sta amr a1 a2 a3 a4 a5]=textread([dir suffix 'ds_' subfault '/_output/fort.gauge'],'%f%f%f%f%f%f%f');
    staG=unique(sta);
    staG=setxor(staG,1016); %Ignore BY.MIC
    Gsta=[];
    for k2=1:Nsta %Loop over all stations for that subfault
        i=find(sta==staG(k2));
        t=a1(i);
        eta=(a5(i));
        tG=0:dt:lt;
        etaG=interp1(t,eta,tG);
        %Now equalize to what teh data for that station actually is
        j=find(numg==staG(k2)); %What is the wave gauge code?
        [tdata etadata]=textread([datapath nameg{j} sta_suffix '.txt'],'%f%f'); %load actual data
        tdata=tdata(2:end); %First sample is inversion weight
        etaG=interp1(tG,etaG,tdata); %Reinterpolate to what is actually int he data
        Gsta=[Gsta ; etaG]; %Append current station
    end
    Gwv(:,2*k)=Gsta; %Add to main matrix
end

%Write to file
save(savef,'Gwv','tGF')
