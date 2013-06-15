function varargout=rtokadawv(workpath,outpath,runName,runID,stasuffix,lambda,G,GSF,Gwv,stations,t_gauges,gauges,dataflag)

% 05/2013 (D.Melgar)
%
% A static slip inversion routine with laplacian regularization and bounds on the lateral and bottom edges of a fault.
% This is a heavily edited version of B.Crowell's rtokada routine modified to include wave gauges in the inversion process 
% and to compute smoothing metrics like Akaike's information criterion and the generalized cross-validation parameter. Also,
% it has been generally  streamlined and optimized for batch runs. THIS DEPRECATES THE ORIGINAL rtokada() FUNCTION.
%
% INPUT VARS
% workpath - Where are teh necessary files for the run held at
% outpath - Where will the output be written to
% runName - String identifying the run name
% runID - Number of this particular run
% stasuffix - Suffix to append to station filename at time of read
% lambda - Smoothing parameter
% G - Green functions for coseismic offsets
% GSF - Green fuctions for seaflorr displacements
% Gwv - Green functions for wave gauges
% stations - File containing coordinates of coseismic offset stations
% t_gauges - Time vector for wave gauges
% gauges - Vector containing name of wave gauge stations
% dataflag - [usegps usewave] determines which data is sued in inversion
%
% OUTPUT VARS
% varargout{1}=lambda - Smoothing parameter used
% varargout{2}=L2 - Inversion misfit || Gm - d ||
% varargout{3}=LS - Solution semi-norm || Lm ||
% varargout{4}=Mo - Solution moment
% varargout{5}=Mw - Solution moment magnitude
% varargout{6}=GCV - Generalized cross-validation parameter
% varargout{7}=AIC - Akaike's information criterion 
% varargout{8}=VRgps - Variance reduction on coseismic offsets
% varargout{9}=RMSwv - RMS misfit onw ave gauges



%SETUP
format long
cd(workpath)
usegps=dataflag(1);  %Which data to fit
usewave=dataflag(2);


%LOAD FAULT MODEL
f2=load('faults_def_small.txt'); %No of fault elements
ast=f2(1);%along strike elements
adi=f2(2);%along dip elements
[xs,ys,zs,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=...
    textread('small_fault.dat','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
[site latinv,loninv]=textread(stations,'%f %f %f');
[latsf lonsf]=textread('seafloor.xy','%f%f');


%REGULARIZATION
%This section computes the regularization matrix to be appended to the bottom of the Green's Function matrix.  It aims to 
%reduce the laplacian between individual fault segments.  To find the nearest segments, it computes the total distance between the
%center of each fault segment and finds the ones "touching" each other
NW=adi-1;
NL=ast-1; %No of fault patches
k=1;
T=[];
for j = 1:NW
    for i = 1:NL
        for m = 1:2
            index1 = (j-1)*NL+i;
            index2 = (j-1)*NL+i-1;
            index3 = (j-1)*NL+i+1;
            index4 = (j-2)*NL+i;
            index5 = (j)*NL+i;
            dx = max(len)/len(index1);
            dy = max(width)/width(index1);
            if (index1 >= 1 && index1 <= length(xs))
                T(k,2*(index1-1)+m) = -2*(dx^-2+dy^-2);
            end
            if (index2 >= 1 && index2 <= length(xs))
                T(k,2*(index2-1)+m) = dx^-2;
            end
            if (index3 >= 1 && index3 <= length(xs))
                T(k,2*(index3-1)+m) = dx^-2;
            end
            if (index4 >= 1 && index4 <= length(xs) )
                T(k,2*(index4-1)+m) = dy^-2;
            end
            if (index5 >= 1 && index5 <= length(xs))
                T(k,2*(index5-1)+m) = dy^-2;
            end
            k=k+1;
        end
    end
end
[h1,h2]=size(T);
Tzeros=zeros(h1,1);

%CREAT MAX SLIP BOUNDS
%Create lower and upper bounds for each fault segment.  It is set up to lock the sides and bottom, and allows the rest of the fault 
%segments to move up to 100 meters in either direction. For positivity, the lower or upper bound needs to be set to zero,
%depending on the sign of slip/strike
eb=1e-2;%maximum displacement the sides and bottom can move, in m
kk=1;
for j=1:NW
    for k=1:NL
        if (k == 1 || k == NL || j == NW)
            lb(kk,1)=-eb;
            ub(kk,1)=eb;
            lb(kk+1,1)=-eb;
            ub(kk+1,1)=eb;
        else
            %These are for strike slip
            lb(kk,1)=-100;
            ub(kk,1)=100;
            %These are for dip slip (mostly thrust)
            lb(kk+1,1)=-1;
            ub(kk+1,1)=100;
        end
        kk=kk+2;
    end
end

%READ DATA
Ugps=[]; %Data vector
Sxyz_gps=[];
%Read land station coseismic offsets
for j = 1:length(site)
    siter = site(j);
    siterb=num2str(siter);
    %y is north, x is east, z is up
    [t,dy,dx,dz,sy,sx,sz]=textread(['neufiles/' siterb '.' stasuffix '.txt'],'%f %f %f %f %f %f %f');
    ux=dx;
    uy=dy;
    uz=dz;
    %Normalize the data weights
    sd=min([sx sy sz]);
    sx=(sx/sd);
    sy=(sy/sd);
    sz=(sz/sd);
    Sxyz_gps=[Sxyz_gps ; sx ; sy ; sz];
    Ucurrent = [ux;uy;uz];
    Ugps=[Ugps;Ucurrent];
end

%Read wave gauges
Ueta=[];
Sxyz_eta=[];
for j=1:length(gauges)
    %First entry is weight, rest is time series
    [tg eta]=textread(['gauges/' gauges{j} '.txt'],'%f%f');
    eta=eta(2:end);
    w=(ones(size(eta))/tg(1));
    Sxyz_eta=[Sxyz_eta ; w];
    Ueta=[Ueta ; eta];
end

% INVERSION
Sxyz=[];
if usegps==1   %Use coseismic offsets
    Ginv = G; 
    U=Ugps;
    Sxyz=[Sxyz ; Sxyz_gps];
else
    Ginv=[];
    U=[];
end
if usewave==1   %Use wave gauge measurements
    Ginv=[Ginv ; Gwv];
    U=[U ; Ueta];
    Sxyz=[Sxyz ; Sxyz_eta];
end  
W=diag(1./Sxyz);
Uinv=[W*U;Tzeros];   %Add weights and zeros for regularization
Gforward=Ginv; %Won't me modified further
Gall=[G ; Gwv]; %GFs used to compute the fits to the different data sets even if they aren't used in the inversion
%Apply weights to data and add zero rows

Tinv=T*lambda; %Apply smoothing parameter
Ginv=[W*Ginv;Tinv];%append the regularization onto the greens function
S=lsqlin(Ginv,Uinv,[],[],[],[],lb,ub);%solve for fault motions, in mm
Uforward = Gforward*S;%%%forward model with original green's function
ngps=size(G,1);
%Get post-inversion metrics
%get L2 norm of misfit
L2=norm(W*U-W*Uforward,2);
%get seminorm of Model
LS=norm(T*S,2);
%get generalized cross validation value
ndata=length(U);
GW=W*Gforward;
Gsharp=(GW'*GW+Tinv'*Tinv)\GW';
GCV=(ndata*(L2^2))/(trace(eye(size(GW*Gsharp))-GW*Gsharp)^2);
%Get Akaike
Ms=max(size(T));
N=max(size(S));
Nhp=2;
phi=(ndata+Ms-N);
AIC=phi*log10(2*pi)+phi*log10(L2^2+(lambda^2)*(LS^2))-2*Ms*log10(lambda)+2*phi*log10(phi)+log10(norm(GW'*GW+(lambda^2)*(T'*T),2))+phi+2*Nhp;
%Now split into GPS and wave gauge metrics
Uforward=Gall*S;
Uforward_gps=Uforward(1:ngps);
Uforward_wv=Uforward(ngps+1:end);
VRgps=sum((Ugps-Uforward_gps).^2)/sum(Ugps.^2);
VRgps=(1-VRgps)*100;
RMSwv=(sum((Uforward_wv-Ueta).^2)/length(Uforward_wv)).^0.5;
%Sea floor displacements
USF=GSF*S;
%Split into strike-slip and dip-slip to compute rake etc.
for i = 1:length(xs)
    S1(i,1)=S((i-1)*2+1)*1000;
    S2(i,1)=S((i-1)*2+2)*1000;
end
ST = (S1.^2+S2.^2).^0.5./1000; %Strike
Mo = sum(30e9.*ST.*area.*1000.*1000)/1e-7;
Mw=(2/3)*log10(Mo)-10.7;
Mo=Mo/1e7;
%And rake
sigma=-(strike-180); %Angle the dipslip vector points to (clockwise is negative) with respect to horizontal
theta=rad2deg(-atan(S1./S2)); %Angle teh slip vector makes in the ss,ds reference frame (clockwise is negative)
rake_h=sigma+theta; %Angle the slip vector make with respect tot eh horizontal due East (clockwise is engative)
fact=5;
rake_amp=ones(size(rake_h))/fact;%sqrt(S1.^2+S2.^2)/fact;
%Output variables
varargout{1}=lambda;
varargout{2}=L2;
varargout{3}=LS;
varargout{4}=Mo;
varargout{5}=Mw;
varargout{6}=GCV;
varargout{7}=AIC;
varargout{8}=VRgps;
varargout{9}=RMSwv;

%Write model results
fid=fopen([outpath runName '.' runID '.slip'],'wt');
fid2=fopen([outpath runName '.' runID '.disp'],'wt');
fid3=fopen([outpath runName '.' runID '.sflr'],'wt');
fid4=fopen([outpath runName '.' runID '.dtopo'],'wt');
fid5=fopen([outpath runName '.' runID '.rake'],'wt');
fid6=fopen([outpath runName '.' runID '.wave'],'wt');
for i=1:length(xs)
    fprintf(fid,'%1.5f %1.5f %1.2f %1.5f %1.5f %1.5f\n',S1(i),S2(i),Mw,xs(i),ys(i),zs(i));
end
nnn = 1;
%Write displacements, observed and synthetic
for i=1:length(site)
    siter = site(i);
    siterb=num2str(siter);
    fprintf(fid2,'%s %1.0f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n',siterb,k,latinv(i),loninv(i),Ugps((nnn-1)*3+2),...
        Ugps((nnn-1)*3+1),Ugps((nnn-1)*3+3),Uforward((nnn-1)*3+2),Uforward((nnn-1)*3+1),Uforward((nnn-1)*3+3));
    nnn = nnn+1;
end
%Write  3 component synthetic seafloor displacements
nnn=1;
for i=1:length(lonsf)
    siter = i;
    siterb=num2str(siter);
    fprintf(fid3,'%s %1.0f %1.4f %1.4f %1.4f %1.4f %1.4f\n',siterb,k,latsf(i),lonsf(i),USF((nnn-1)*3+2),USF((nnn-1)*3+1),USF((nnn-1)*3+3));
    nnn = nnn+1;
end
%Write  GeoClaw dtopo type 3 file
for i=1:length(lonsf)
    fprintf(fid4,'%1.0f %1.4f %1.4f %1.4f\n',0,lonsf(i),latsf(i),0);
end
nnn=1;
for i=1:length(lonsf)
    fprintf(fid4,'%1.0f %1.4f %1.4f %1.4f\n',1,lonsf(i),latsf(i),USF((nnn-1)*3+3));
    nnn = nnn+1;
end
%Write rake vector information
for i=1:length(xs)
    
    fprintf(fid5,'%1.4f %1.4f %1.4f %1.4f\n',xs(i),ys(i),rake_h(i),rake_amp(i));
end
%Write wave gauges, observed and synthetic: time,observed,synthetic
nnn=1;
for i=1:length(t_gauges)
    fprintf(fid6,'%1.4f %1.4f %1.4f\n',t_gauges(i),Ueta(i),Uforward(ngps+i));
end
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);

%Output to screen
display(['   lambda = ' num2str(lambda)])
display(['   VRgps = ' num2str(VRgps) '%'])
display(['   RMSwv = ' num2str(RMSwv) ''])
display(['   || Gm-d || = ' num2str(L2)])
display(['   || LM || = ' num2str(LS)])
display(['   GCV = ' num2str(GCV)])
display(['   AIC = ' num2str(AIC)])
display(['   Mw = ' num2str(Mw)])