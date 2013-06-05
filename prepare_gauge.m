function prepare_gauge

%DMM 05/2013
%
%Prepare wave gauges for ingestion into inverse problem

maxt=50*60;
dt=15;
outdir='/Users/dmelgarm/Research/Data/Tohoku/RTOkada/gauges/';
obsdir='/Users/dmelgarm/Research/Data/Tohoku/Gauges';
suffix='_50min';
gaugedir='/Users/dmelgarm/Research/Data/Tohoku/Gauges';
gauges=[1002 1004 1005 1009 1010 1011 1015 1016 1017 1018 1019 1020 1021 1023 1028 1029];
weights=[1 1 5 5 5 5 1 5 5 1 1 1 1 1 10 10];
%Earthquake origin time
t0=11*86400+5*3600+46*60;
%dt in synthetics
dtsynth=15; %secs
%Station list and correspondences
sta_obs={'TI.AYU','TI.CHO','BY.FUK','BY.IWC','BY.IWN','BY.IWS','TI.MER','BY.MIC','BY.MIN','TI.MYE','TI.MYO','TI.OFU','TI.OKA','TI.ONA','OB.TM1','OB.TM2'};
ti=0:dt:maxt;
for k=1:length(sta_obs)
    sta_obs{k}
    %Get observed
    cd(obsdir)
    [day a a hr min sec etao]=textread([sta_obs{k} '.txt'],'%f%f%f%f%f%f%f');
    to=day*86400+hr*3600+min*60+sec;
    to=to-t0;
    etai=interp1(to,etao,ti);
    %Clean up nan's
    i=find(~isnan(etai));
    etaout=[weights(k) etai(i)-etai(i(1))];
    tout=[weights(k) ti(i)];
    if k==8 %Manually edit BY.IWC
        i=find(tout<=1485);
        tout=tout(i);
        etaout=etaout(i);
    end
    plot(tout(2:end),etaout(2:end))
    title(sta_obs{k})
    pause()
    out=[tout' etaout'];
    save([outdir sta_obs{k} suffix '.txt'],'out','-ascii')
end