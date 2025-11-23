% created by Chen Xintong on 2025-11-23

clear;
clc;
close all;

addpath("..");
addpath("../Geometry");
addpath("../ERP");
addpath("../Pathloss");

para_init;

%% load data
load('../../data/ricianProve/2025-11-23_channelPrim.mat');
path = [3 30];

%% plot
% with ERP
linewidth_theory = 6;
linewidth_MC = 5;
MarkerSize_MC = 10;
x_temp{1} = 10*log10(x_temp{1}.^2);
x_temp{2} = 10*log10(x_temp{2}.^2);
index = [100:300:2*length(Channelprim_MC(:,1,1))+1 5999];
index2 = [10:300:2*length(Channelprim_MC(:,1,1))+1 5999];

scheme1 = 10*log10(mean(abs(Channelprim_MC(:,2,1)).^2));
scheme1 = eval(vpa(scheme1,3));
scheme2 = 10*log10(mean(abs(Channelprim_MC(:,2,2)).^2));
scheme2 = eval(vpa(scheme2,3));

figure(1);
set(gcf,'color','w')
h1 = plot(x_temp{1},Channelprim_Theory{1});
set(h1,'linestyle','-','linewidth',linewidth_theory,'color',[1.00,0.32,0.41])
hold on
h2 = cdfplot(10*log10(abs(Channelprim_MC(:,1,1)).^2));
h2.XData = h2.XData(index);
h2.YData = h2.YData(index);
set(h2,'linestyle',':','linewidth',linewidth_MC,'Marker','^','MarkerSize',MarkerSize_MC,'MarkerFaceColor',color(6,:),'color',color(6,:));
hold on
h3 = cdfplot(10*log10(abs(Channelprim_MC(:,2,1)).^2));
h3.XData = h3.XData(index);
h3.YData = h3.YData(index);
set(h3,'linestyle','-.','linewidth',linewidth_MC,'Marker','s','MarkerSize',MarkerSize_MC,'MarkerFaceColor','k','color','k');
hold on

hline1 = plot([scheme1,scheme1],[0,0.57],'linestyle','--','linewidth',2,'color','k');
text(scheme1, 0, num2str(scheme1),'fontsize',28);
set(get(get(hline1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h4 = plot(x_temp{2},Channelprim_Theory{2});
set(h4,'linestyle','-','linewidth',linewidth_theory,'color',color(3,:));
hold on
h5 = cdfplot(10*log10(abs(Channelprim_MC(:,1,2)).^2));
h5.XData = h5.XData(index2);
h5.YData = h5.YData(index2);
set(h5,'linestyle',':','linewidth',linewidth_MC,'Marker','^','MarkerSize',MarkerSize_MC,'MarkerFaceColor','g','color','g');
hold on
h6 = cdfplot(10*log10(abs(Channelprim_MC(:,2,2)).^2));
h6.XData = h6.XData(index2);
h6.YData = h6.YData(index2);
set(h6,'linestyle','-.','linewidth',linewidth_MC,'Marker','s','MarkerSize',MarkerSize_MC,'MarkerFaceColor',[0.89,0.35,1.00],'color',[0.89,0.35,1.00]);

hline2 = plot([scheme2,scheme2],[0,0.52],'linestyle','--','linewidth',2,'color','k');
text(scheme2, 0,num2str(scheme2),'fontsize',28);
set(get(get(hline2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% hold on
% plot([-80,x_temp{1}(ZeroPointFiveIdx),x_temp{2}(ZeroPointFiveIdx1)],[0.5 0.5 0.5],'linestyle','--','linewidth',2,'color','k')

axis([-80,-55,0,1])
title('')
xlabel("channel power gain (dB)","interpreter","latex")
ylabel("CDF of channel power gain $\left| h^{\prime}_{\rm{I}_{0} 1} \right|^2$","interpreter","latex")
set(gca,'fontsize',32)
MCcurve1 = ['$\theta_\textrm{g}=\pi/12$',', $n_\textrm{p}$=',num2str(path(1)),', MC'];
MCcurve2 = ['$\theta_\textrm{g}=\pi/12$',', $n_\textrm{p}$=',num2str(path(2)),', MC'];
MCcurve3 = ['$\theta_\textrm{g}=5\pi/12$',', $n_\textrm{p}$=',num2str(path(1)),', MC'];
MCcurve4 = ['$\theta_\textrm{g}=5\pi/12$',', $n_\textrm{p}$=',num2str(path(2)),', MC'];
legend("$\theta_\textrm{g}=\pi/12$, Formula",MCcurve1,MCcurve2,"$\theta_\textrm{g}=5\pi/12$, Formula",MCcurve3,MCcurve4,'Location','best','interpreter','latex','fontsize',35)

%%% ZOOM
axes('Position',[0.74 0.37 0.25 0.25]);
box on;

h1 = plot(x_temp{1},Channelprim_Theory{1});
set(h1,'linestyle','-','linewidth',linewidth_theory,'color',[1.00,0.32,0.41])
hold on
h2 = cdfplot(10*log10(abs(Channelprim_MC(:,1,1)).^2));
h2.XData = h2.XData(index);
h2.YData = h2.YData(index);
set(h2,'linestyle',':','linewidth',linewidth_MC,'Marker','^','MarkerSize',MarkerSize_MC,'MarkerFaceColor',color(6,:),'color',color(6,:));
hold on
h3 = cdfplot(10*log10(abs(Channelprim_MC(:,2,1)).^2));
h3.XData = h3.XData(index);
h3.YData = h3.YData(index);
set(h3,'linestyle','-.','linewidth',linewidth_MC,'Marker','s','MarkerSize',MarkerSize_MC,'MarkerFaceColor','k','color','k');

xlim([-66, -65]);
ylim([0.6, 0.7]);
title('');
xlabel('');
ylabel('');
set(gca,'fontsize',24);
