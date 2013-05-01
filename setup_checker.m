function setup_checker

workpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada'
suffix='1patch'
suffix='checker_kal_gps';
%Make checkerboard
slip=10; %in m
%Ch=zeros(1600,1);
Ch=zeros(378,1);

% %Small checkerboard
% S=2:4:78;
% j=S;
% for k=1:19
%     if mod(k,2)==0
%         S=[S (j+80*(k)-1)];
%     else
%         S=[S (j+80*(k)+1)];
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Big checkerboard
% S1=2:8:158;
% S2=4:8:158;
% S=interleave(S1,S2);
% clear S1 S2
% j=S;
% for k=1:9
%     if mod(k,2)==0
%         S=[S (j+160*(k))];
%     else
%         S=[S (j+160*(k)+4)];
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Single fault patch (Big)
% S=280;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single fault patch (small)
%S=150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Three asperity checekrboard (Big fault)
% S1=22:2:60;
% S2=S1+80;
% S3=S2+80;
% S=[S1 S2 S3];
% clear S2 S3
% S1=S1+7*80;
% S2=S1+80;
% S3=S2+80;
% S=[S S1 S2 S3];
% clear S2 S3
% % S1=S1+7*80;
% % S2=S1+80;
% % S3=S2+80;
% % S=[S S1 S2 S3];
% clear S1 S2 S3;

%Two asperity checekrboard (small fault)
S1=10:2:34;
S2=S1+42;
S=[S1 S2];
clear S2
S1=S1+210;
S2=S1+42;
S=[S S1 S2];
clear S1 S2;


Ch(S)=slip;
%Load GFs
cd(workpath)
%load green.mat
load green_small_kal.mat
%Make displacements
D=G*Ch;
%Wrie to file
%Load all GPS
% cd ..
% load tohoku_gps_all
% coseis2=coseis;
% load  tohoku_gps_statics
%Load Kalman
%Process file sfor 138 kalman stations
cd ..
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
%

cd(workpath)
N=length(coseis.lon);
cd neufiles
kt=0;
for k=1:N
    if (isnan(coseis.Nkal(k))==1 || isnan(coseis.Ekal(k))==1 || isnan(coseis.Ukal(k))==1)
        a=0;
    else
        kt=kt+1;
        coseis.gname(k)
        P=[1 D(3*kt-1) D(3*kt-2) D(3*kt) coseis2.stdn(k) coseis2.stde(k) coseis2.stdu(k)];
        save([num2str(coseis.gname(k)) '.' suffix '.txt'],'P','-ascii')
    end
end
cd ..
%Now wriote checkerboard in RTO format
[xs,ys,zs,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=textread('small_fault.dat','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
for i = 1:length(xs)
    S1(i,1)=Ch((i-1)*2+1)*1000;
    S2(i,1)=Ch((i-1)*2+2)*1000;
end
%DEBUG
lon=138:0.1:148;
lat=33:0.1:44;
lat=fliplr(lat);
k=0;
for i=1:length(lat)
    for j=1:length(lon)
        k=k+1;
        SF(k,1)=lat(i);
        SF(k,2)=lon(j);
    end
end
% END DEBUG
ST = (S1.^2+S2.^2).^0.5./1000;
Mo = sum(30e9.*ST.*area.*1000.*1000)/1e-7;
Mw=2/3*log10(Mo)-10.7
Mo=Mo/1e7;
if (Mo == 0)
    Mw = 0;
end
cd output
fid=fopen('checkerboard_2patch.0001.slip','wt');
for i=1:length(xs)
    fprintf(fid,'%1.0f %1.0f %1.5f %1.5f %1.2f %1.5f %1.5f %1.5f %1.4f\n',1,i,S1(i),S2(i),Mw,xs(i),ys(i),zs(i),0);
end
cd(workpath)
fclose(fid)
a=0;


