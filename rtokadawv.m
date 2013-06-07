function varargout=rtokadawv(workpath,outpath,runName,runID,stasuffix,lambda,G,GSF,Gwv,weightflag,stations,t_gauges,gauges)

% 05/2013 (DM)
%
% This is a modified version fo thr rtokada() routine to incldue support for inversion of wavegauges as well



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Major Variables and Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
cd(workpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create the fault planes, from faults.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f2=load('faults_defl.txt'); %No of fault elements
f2=load('faults_def_small.txt'); %No of fault elements
ast=f2(1);%along strike elements
adi=f2(2);%along dip elements
[xs,ys,zs,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=textread('small_fault.dat','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
[site latinv,loninv]=textread(stations,'%f %f %f');
[latsf lonsf]=textread('seafloor.xy','%f%f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Greens Function Calculation, faults file is as follows
%Coordiantes of the centre of the patch, coordinates of the 4
%cornerns,strike,dip,along strike elngth,along dip length, total area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[xs,ys,zs,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,strike,dip,len,width,area]=textread('fault_plane.dat','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

%
% tic
% [G] = okada_green(xrs,yrs,zrs,strike,dip,width.*1000,len.*1000);
% toc
% save('green_small_kal.mat','G');


% tic
% [GSF] = okada_green(xrs,yrs,zrs,strike,dip,width.*1000,len.*1000);
% toc
% save('greenSF_small_kal.mat','GSF');
% a=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section computes the regularization matrix to be appended to the
%bottom of the Green's Function matrix.  It aims to reduce the laplacian
%between individual fault segments.  To find the nearest segments, it
%computes the total distance between the center of each fault segment
%and finds the ones "touching" each other
NW=adi-1;
NL=ast-1; %No of fault patches
k=1;
T=[];
for j = 1:NW%put 1 to clamp top or 2 to unclamp top
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
TU=zeros(h1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create lower and upper bounds for each fault segment.  It is set up to
%lock the sides and bottom, and allows the rest of the fault segments to
%move up to 10 meters in either direction
%For positivity, the lower or upper bound needs to be set to zero,
%depending on the sign of slip/strike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            %             lb(kk+1,1)=-1;  %For checkerboard stuff
            %             ub(kk+1,1)=2000;
        end
        kk=kk+2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write data tot he following files
fid=fopen([outpath runName '.' runID '.slip'],'wt');
fid2=fopen([outpath runName '.' runID '.disp'],'wt');
fid3=fopen([outpath runName '.' runID '.sflr'],'wt');
fid4=fopen([outpath runName '.' runID '.dtopo'],'wt');
fid5=fopen([outpath runName '.' runID '.rake'],'wt');
fid6=fopen([outpath runName '.' runID '.wave'],'wt');


% ydist = (loninv-lono)*llon/1000;
% xdist = (latinv-lato)*llat/1000;
% tdist = (xdist.^2+ydist.^2).^0.5;

%%%%%%
%This is the weight matrix, which is currently not being used
k=1;
% for i = 1:length(site)
%     C(k,k)=min(tdist)^2/tdist(i)^2/10;
%     C(k+1,k+1)=min(tdist)^2/tdist(i)^2/10;
%     C(k+2,k+2)=min(tdist)^2/tdist(i)^2/50;
%     k=k+3;
% end



GS = G;
kk=1;
for k = 1:1
    G2 = G;
    GS2 = GS;
    %C2 = C;
    UD=[];
    Sxyz=[];
    StaRem=[];
    for j = 1:length(site)
        siter = site(j);
        %sitera = char(siter);
        %siterb = sprintf('%s',sitera);
        siterb=num2str(siter);
        %y is north, x is east, z is up
        [t,dy,dx,dz,sy,sx,sz]=textread(['neufiles/' siterb '.' stasuffix '.txt'],'%f %f %f %f %f %f %f');
        a1 = find (k == t);
        if (~isempty(a1))
            %             uxx = movavg(dx,ma*100,ma*100);%100 is for the rate (100 Hz)
            %             uyy = movavg(dy,ma*100,ma*100);
            %             uzz = movavg(dz,ma*100,ma*100);
            %             ux = uxx(a1);
            %             uy = uyy(a1);
            %             uz = uzz(a1);
            ux=dx;
            uy=dy;
            uz=dz;
            Sxyz=[Sxyz ; sx ; sy ; sz];
            U1 = [ux;uy;uz];
            UD = [UD;U1];
        else
            StaRem = [StaRem;j];
        end
    end
    %Normalize weight of land GPS stations to a 1 %Maybe try weighting by
    %the eman of the noise??
    gpsmult=1;
    Stemp=1./Sxyz;
    Stemp=Stemp./max(Stemp);
    Stemp=Stemp*gpsmult;
    Sxyz=1./Stemp;
    %Read gauges
    for j=1:length(gauges)
        %First entry is weight, rest is time series
        [tg eta]=textread(['gauges/' gauges{j} '.txt'],'%f%f');
        eta=eta(2:end);
        w=ones(size(eta))/tg(1);
        Sxyz=[Sxyz ; w];
        UD=[UD ; eta];
    end
    %Add tsunami GFs
    G2=[G2 ; Gwv];
    GS2=[GS2 ; Gwv];
    %Weight by station noise
    if weightflag==1
        W=diag(1./Sxyz);
    else
        W=eye(length(site)*3);
    end
    %%%%%%
    
    if (~isempty(UD))
        UD2=[W*UD;TU];
        
        T2 = T*lambda;%Place smoothness constraint here
        G3=[W*G2;T2];%append the regularization onto the greens function
        tic
        S=lsqlin(G3,UD2,[],[],[],[],lb,ub);%solve for displacements, in mm
        toc
        %S=lsqlin(G,UD,[],[],[],[],[],[]);%solve for displacements, in mm
        %S=G\UD;
        %         [U s V]=cgsvd(G,T);   %Get SVD
        %         [S,rho,eta]=tikhonov(U,s,V,UD,lambda);
        %         S0=ones(length(lb),1);
        %         S = l1decode_pd(S0, G3, [], UD2);
        UP = GS2*S;%%%forward model with original green's function
        ngps=size(G,1)/3;
        %Divide into GPS and wave gauge data
        UPgps=UP(1:ngps);
        UDgps=UD(1:ngps);
        UPwv=UP(ngps+1:end);
        UDwv=UD(ngps+1:end);
        %Get post-inversion metrics
        rms=(sum((UP-UD).^2)/length(UP)).^0.5;
        %get L2 norm of misfit
        L2=norm(W*UP-W*UD,2);
        %get seminorm of Model
        LS=norm(T*S,2);
        %get generalized cross validation value
        ndata=length(UD);
        GW=W*GS2;
        Gsharp=(GW'*GW+T2'*T2)\GW';
        GCV=(ndata*(L2^2))/(trace(eye(size(GW*Gsharp))-GW*Gsharp)^2);
        %Get Akaike
        Ms=max(size(T2));
        N=max(size(S));
        Nhp=2;
        phi=(ndata+Ms-N);
        ABIC=phi*log10(2*pi)+phi*log10(L2^2+LS^2)-2*Ms*log10(lambda)+2*phi*log10(phi)+log10(norm(GW'*GW+T2'*T2,2))+phi+2*Nhp;
        %Now split into GPS and wave gauge metrics
        VRgps=sum((UDgps-UPgps).^2)/sum(UDgps.^2);
        VRgps=(1-VRgps)*100;
        VRwv=sum((UDwv-UPwv).^2)/sum(UDwv.^2);
        VRwv=(1-VRwv)*100;
        %Sea floor displacements
        USF=GSF*S;
        
        for i = 1:length(xs)
            S1(i,1)=S((i-1)*2+1)*1000;
            S2(i,1)=S((i-1)*2+2)*1000;
        end
        ST = (S1.^2+S2.^2).^0.5./1000;
        Mo = sum(30e9.*ST.*area.*1000.*1000)/1e-7;
        Mw=(2/3)*log10(Mo)-10.7;
        Mo=Mo/1e7;
        VR=sum((UD-UP).^2)/sum(UD.^2);
        VR=(1-VR)*100;
        if (Mo == 0)
            Mw = 0;
        end
        %And rake
        sigma=-(strike-180); %Angle the dipslip vector points to (clockwise is negative) with respect to horizontal
        theta=rad2deg(-atan(S1./S2)); %Angle teh slip vector makes in the ss,ds reference frame (clockwise is negative)
        rake_h=sigma+theta; %Angle the slip vector make with respect tot eh horizontal due East (clockwise is engative)
        fact=5;
        rake_amp=ones(size(rake_h))/fact;%sqrt(S1.^2+S2.^2)/fact;
        %Output
        varargout{1}=lambda;
        varargout{2}=L2;
        varargout{3}=LS;
        varargout{4}=Mo;
        varargout{5}=Mw;
        varargout{6}=VR;
        varargout{7}=GCV;
        varargout{8}=ABIC;
        varargout{9}=VRgps;
        varargout{10}=VRwv;
        
        %Write model results
        for i=1:length(xs)
            fprintf(fid,'%1.0f %1.0f %1.5f %1.5f %1.2f %1.5f %1.5f %1.5f %1.4f\n',k,i,S1(i),S2(i),Mw,xs(i),ys(i),zs(i),rms);
        end
        nnn = 1;
        %Write displacements, observed and synthetic
        for i=1:length(site)
            siter = site(i);
            siterb=num2str(siter);
            fprintf(fid2,'%s %1.0f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n',siterb,k,latinv(i),loninv(i),UD((nnn-1)*3+2),UD((nnn-1)*3+1),UD((nnn-1)*3+3),UP((nnn-1)*3+2),UP((nnn-1)*3+1),UP((nnn-1)*3+3));
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
        %Write wave gauges, observed and synthetic time,observed,synthetic
        nnn=1;
        for i=1:length(t_gauges)
            fprintf(fid6,'%1.4f %1.4f %1.4f\n',t_gauges(i),UD((3*ngps)+i),UP((3*ngps)+i));
        end
    end
end
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
%Output to screen
display(['   lambda = ' num2str(lambda)])
display(['   VR = ' num2str(VR) '%'])
display(['   VRgps = ' num2str(VRgps) '%'])
display(['   VRwv = ' num2str(VRwv) '%'])
display(['   || LM || = ' num2str(LS)])
display(['   GCV = ' num2str(GCV)])
display(['   ABIC = ' num2str(ABIC)])
display(['   Mw = ' num2str(Mw)])