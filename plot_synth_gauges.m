function plot_synth_gauges(num)

%DMM 06/2013
%
% Plot and compare wave gauge results of the inversion

%read data
runName='RTOkada_kalwv10min';
outpath='/Users/dmelgarm/Research/Data/Tohoku/RTOkada/output/';
[t obs syn]=textread([outpath runName '.' num '.wave'],'%f%f%f');
cd(outpath)
load('../tohoku_wvGF_60min.mat')%This is needed to name the stations
%Set up base plot
figure
set(gcf,'Position',[1921 -41 1600 1127]);
ha = tight_subplot(4, 4, [0 0.05], 0.1, 0.1);
tlim=3600;
%Parse data vectors and plot
for k=1:length(gauges)
    i=find(diff(t)<15);
    i=[1 ; (i+1) ; length(t)];
    %Get observed
    to=t(i(k):i(k+1)-1);
    etaobs=obs(i(k):i(k+1)-1);
    %Get synthetic
    ts=to;
    etasyn=syn(i(k):i(k+1)-1);
    axes(ha(k))
    plot(ts/60,etaobs,ts/60,etasyn,'LineWidth',2);
    %plot(ts/60,etas)
    set(gca,'LineWidth',1.5);
    xlim([0 tlim/60]);
    yl=ylim;
    staname=[gauges{k}];
    text(10,0.8*yl(2),staname(1:6));
    if k==1
        legend('Observed','Modeled')
    end
    if (k==1 | k ==5 | k==9 | k==13 | k==17)
        ylabel('\eta(m)')
    end
    if (k==13 | k==14 | k==15 | k==16)
        xlabel('Minutes after OT');
    end
    if (k < 13)
        set(gca,'XTickLabel',[]');
    end
end