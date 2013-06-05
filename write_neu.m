function write_neu

%DMM 03/2013
%
% Read in data from standard matlab structure and output to NEU format

%suffix='raw';
suffix='kal_weightgps'
cd('/Users/dmelgarm/Research/Data/Tohoku')
%load tohoku_gps_all
%coseis2=coseis;
%load  tohoku_gps_statics

%Process file sfor 138 kalman stations
load tohokukal
tcomp=157;
i=find(tohoku.process==1)
for k=1:length(i);
%     %For Kalman data
    itkal=find(tohoku.tmakal(:,i(k))>tcomp-0.001 & tohoku.tmakal(:,i(k))<tcomp+0.001);
    coseis.Nkal(k)=tohoku.nmakal(itkal,i(k));
    coseis.Ekal(k)=tohoku.emakal(itkal,i(k));
    coseis.Ukal(k)=tohoku.umakal(itkal,i(k));
    %For GPS data
    it=find(tohoku.tmagps(:,i(k))>tcomp-0.001 & tohoku.tmagps(:,i(k))<tcomp+0.001);
    coseis.N(k)=tohoku.nmagps(it,i(k));
    coseis.E(k)=tohoku.emagps(it,i(k));
    coseis.U(k)=tohoku.umagps(it,i(k));
end
coseis.lon=tohoku.lon(i);
coseis.lat=tohoku.lat(i);
coseis.gname=tohoku.gname(i);
% %For Kalman
% coseis.stdn=tohoku.stdnkal(i);
% coseis.stde=tohoku.stdekal(i);
% coseis.stdu=tohoku.stdukal(i);
%For GPS
coseis.stdn=tohoku.stdngps(i);
coseis.stde=tohoku.stdegps(i);
coseis.stdu=tohoku.stdugps(i);
coseis2=coseis;
% Doen reading Kalmans tuff
N=length(coseis.lon);
cd RTOkada/neufiles
stak=0;
for k=1:N
    if (isnan(coseis.Nkal(k))==1 || isnan(coseis.Ekal(k))==1 || isnan(coseis.Ukal(k))==1)
        a=0;
    else
        stak=stak+1;
        P=[1 coseis.Nkal(k) coseis.Ekal(k) coseis.Ukal(k) coseis2.stdn(k) coseis2.stde(k) coseis2.stdu(k)];
        save([num2str(coseis.gname(k)) '.' suffix '.txt'],'P','-ascii')
        sta(stak)=coseis.gname(k);
        lat(stak)=coseis.lat(k);
        lon(stak)=coseis.lon(k);
    end
end
%Save styations
P=[sta' lat' lon'];
cd ..
%save('stations.xy','P','-ascii')
save('stations_kal.xy','P','-ascii')
%Places where I want seafloor output
clear P
lon=138:2/60:148;
lat=33:2/60:44;
lat=fliplr(lat);
k=0;
for i=1:length(lat)
    for j=1:length(lon)
        k=k+1;
        P(k,1)=lat(i);
        P(k,2)=lon(j);
    end
end
cd ..
save('seafloor.xy','P','-ascii')