function plot_lcurve

N=200;
%maxmod=550;  %for _wv
maxmod=400;  %for _wv_highGPS
%fcut=50;

%runName='RTOkada'
runName='RTOkada_wv'
%runName='Checkerboard'
cd('/Users/dmelgarm/Research/Data/Tohoku/RTOkada/output')

for k=1:N
    if k<10
        runID=['000' num2str(k)];
    elseif k<100
        runID=['00' num2str(k)];
    elseif k<1000
        runID=['0' num2str(k)];
    else
        runID=num2str(k);
    end
    [a n]=textread([runName '.' runID '.log'],'%s%f');
    lambda(k)=n(2);
    Lm(k)=n(4);
    L2(k)=n(3);
    GCV(k)=n(8);
    ABIC(k)=n(9);
end
i=find(Lm<maxmod);
Lm=Lm(i);
lambda=lambda(i);
L2=L2(i);
% clf
% ha = tight_subplot(2, 2, 0.08, 0.1, 0
% 
% axes(ha(1));
% i=find(Lm<12e10);
% % i=find(lambda>=0.1 & lambda <=100);
% %i=1:200;
%  L2=L2(i);
%  Lm=Lm(i);
%  lambda=lambda(i);
%  GCV=GCV(i);
%  ABIC=ABIC(i);
% loglog(L2,Lm,'ok','MarkerSize',3,'MarkerFaceColor','k')
% hold on
% xlabel('||Gm-d||_2')
% ylabel('||Lm||_2')
% posk=1.05;
% for k=1:1:length(i)
%     text(L2(k)*posk,Lm(k)*posk,['\lambda = ' num2str(lambda(k))],'FontSize',16);
%     plot(L2(k),Lm(k),'*r','MarkerSize',12,'LineWidth',1.7)
% end
% [ppx,kappa,reg_c,rho_c,eta_c] = l_corner(L2',Lm',lambda');
% kappa=abs(kappa);
% hold on
% plot(rho_c,eta_c,'*b','MarkerSize',14,'LineWidth',1.7)
% xl=[min(L2) max(L2)];
% xlim(xl);
% ylim([min(Lm) max(Lm)]);
% 
% axes(ha(3));
% semilogx(10.^ppx,abs(kappa),'b','LineWidth',2)
% xlabel('||Gm-d||_2')
% ylabel('Curvature')
% xlim(xl);
% ylim([min(kappa) max(kappa)]);
% 
% axes(ha(2))
% loglog(lambda,GCV,'LineWidth',2)
% ylabel('V_0(\alpha)')
% xlabel('\lambda')
% display(['lambda = ' num2str(reg_c)]);
% xlim([min(lambda) max(lambda)])
% 
% axes(ha(4))
% semilogx(lambda,ABIC,'m','LineWidth',2)
% ylabel('AIC')
% xlabel('\lambda')
% xlim([min(lambda) max(lambda)])

clf
ha = tight_subplot(2, 1, 0, 0.1, 0.3);
%Smooth the curve
% fcut=5;
% L2=L2(2:end);
% Lm=Lm(2:end);
% x=linspace(log10(min(L2)),log10(max(L2)),200);
% y=interp1(log10(L2),log10(Lm),x);
% yf=fbutter(y,x(2)-x(1),4,fcut,'low');
% L2i=10.^x;
% Lmi=10.^yf;
% %Assign lambdas
% lambdai=interp1(Lm,lambda,Lmi);
% %Update everything
% L2=L2i;
% Lm=Lmi;
% lambda=lambdai;

[ppx,kappa,reg_c,rho_c,eta_c] = l_corner(L2',Lm',lambda');
kappa=abs(kappa);
x=linspace(log10(min(L2)),log10(max(L2)),200);
ykappa=interp1(ppx,kappa,x);
i=find(~isnan(ykappa));
x=x(i);
ykappa=ykappa(i);
%yf=fbutter(ykappa,x(2)-x(1),6,fcut,'low');
yf=ykappa;
kappa=yf;
ppx=x;

axes(ha(1));
i=find(Lm<maxmod);
L2=L2(i);
Lm=Lm(i);
lambda=lambda(i);
Lmi=interp1(log10(L2),log10(Lm),ppx); %I need this so I know whick kapaps to throw out int he itneproalted scheme
lambdai=interp1(log10(L2),log10(lambda),ppx);
i=find(Lmi<maxmod);
kappa=kappa(i);
ppx=ppx(i);
i=find(isnan(ppx));
%Find corner in original data
[kmax imax]=max(kappa);
L2corner=ppx(imax);
[junk i]=min(abs(x-L2corner));
reg_c=10^lambdai(i);
rho_c=10^L2corner;
eta_c=10^Lmi(imax);
kappa=[kappa 0];
ppx=[ppx log10(max(L2))];

loglog(L2,Lm,'ok','MarkerSize',3,'MarkerFaceColor','k')
hold on
ylabel('||Lm||_2')
posk=1.05;
for k=10:20:length(lambda)
    text(L2(k)*posk,Lm(k)*posk,['\lambda = ' num2str(lambda(k))],'FontSize',16);
    plot(L2(k),Lm(k),'*r','MarkerSize',12,'LineWidth',1.7)
end
hold on
plot(rho_c,eta_c,'*b','MarkerSize',14,'LineWidth',1.7)
xl=[min(L2) max(L2)];
xlim(xl);
ylim([min(Lm) max(Lm)]);

axes(ha(2));
semilogx(10.^ppx,abs(kappa),'b','LineWidth',2)
xlabel('||Gm-d||_2')
ylabel('\kappa')
xlim(xl);
ylim([min(kappa) max(kappa)]);


display(['lambda = ' num2str(reg_c)]);