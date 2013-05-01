function convert_fault

%Convert from corner coordinates to center coordinates

cd /Users/dmelgarm/Research/Data/Tohoku/RTOkada
load Romano_etal_fault.mat
%read corners
ur=F(:,1:2);
ul=F(:,3:4);
bl=F(:,5:6);
br=F(:,7:8);
zt=F(:,9);
strike=F(:,10);
dip=F(:,11);
%Get along strike and dip lengths
strike_length=25*ones(size(strike));
dip_length=25*ones(size(dip));
area=strike_length.*dip_length;
%Compute middle (lat and lon)
lonc=(ul(:,1)-br(:,1))/2;
lonc=br(:,1)+lonc;
latc=(ul(:,2)-br(:,2))/2;
latc=br(:,2)+latc;
%Compute middle along top left to bottom left edge
lonedge=(ul(:,1)-bl(:,1))/2;
lonedge=bl(:,1)+lonedge;
latedge=(ul(:,2)-bl(:,2))/2;
latedge=bl(:,2)+latedge;
%Figure out hte depths 9 subfaults at a time
i=1:9;
Li=length(i);
dz=zeros(size(i));
for k=1:21
    j=i+(k-1)*Li;
    zc(j)=zt(j);
    dz=diff(zt(j));
    dz=dz/2;
    dz(9)=0;
    zc(j)=zc(j)+dz';
end
%Now fix z-coordinate of bottom patches
i=9:9:189;
dz=zeros(size(zc));
zbottom=dz;
dz(i)=4;
zc=zc+dz;
%Change order to along strike starting at bottom right corner
i=1:9:189;
order=[fliplr(i)];
for k=1:8
    i=i+1;
    order=[order fliplr(i)];
end
%depth of bottom nodes
i=9:9:189;
i=[i i+189 i+(2*189) i+(3*189)];
zbottom=[[zt(2:end) ;zt(end)+8]; [zt(2:end) ;zt(end)+8];zt ;zt ];
    zbottom(i)=zbottom(i-1)+8;
    

%Write to file
F=[lonc(order) latc(order) -zc(order)'*1000 zeros(size(latc,1),12) strike(order) dip(order) strike_length(order) dip_length(order) area(order)];
save('small_fault.dat','F','-ascii')
%Plot
x=[ur(:,1) ; ul(:,1)];
y=[ur(:,2) ; ul(:,2)];
z=[zt ;zt];
scatter3(x,y,-z)
hold on
scatter3(lonc,latc,-zc,'r')
a=0;