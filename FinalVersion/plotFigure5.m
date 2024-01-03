% created by Chen Xintong on 2023-12-17
clc;
clear;
close all;

para_init;

%% initialized
% load data of Fig.5
load('../data/singleCell/20231217avgThroughputinSingleCell');
load('../data/singleCell/20231217barbarSinSingleCell');

% other variables
nnn = 3;
[UE_2D_raw, UE_2D, sel_idx, ~, ~, ~, ~] = calLocationAndAngle(nnn);
numUEs = size(UE_2D,1); % number of UEs

%% plot
% Average Ergodic Throughput (AET) bps/Hz
AETRate(1,:) = sum(avgThroughput(:,:,1)) / numUEs; % Ground UE
AETRate(2,:) = sum(avgThroughput(:,:,2)) / numUEs; % Aerial UE

% Fairness J
Fairness(1,:) = sum(avgThroughput(:,:,1)).^2 ./ ( length(avgThroughput) .* sum(avgThroughput(:,:,1).^2) ); % Ground UE
Fairness(2,:) = sum(avgThroughput(:,:,2)).^2 ./ ( length(avgThroughput) .* sum(avgThroughput(:,:,2).^2) ); % Aerial UE

% Average Mean Signal Power (dB)
AvgMeanSignalPower(1,:) = 10.*log10( mean(barbarS(:,:,1),1) ); % Ground UE
AvgMeanSignalPower(2,:) = 10.*log10( mean(barbarS(:,:,2),1) ); % Aerial UE

AETRate
Fairness
AvgMeanSignalPower

% Ergodic Throughput distribution plot
temp_ue_x = 0:3:2*r; % x
temp_ue_y = -sqrt(3)/2*r:10:sqrt(3)/2*r; % y

%% Proposed schemes
figure(2);
set(gcf,'Color','w')
[Xmesh3D,Ymesh3D,Zmesh3D] = meshgrid(temp_ue_x,temp_ue_y,Hue);
for nis = 1:numIRSscheme
    for np = 1:numPlot
        temp = NaN*ones(length(UE_2D_raw),1);
        temp(sel_idx) = avgThroughput(:,nis,np);
        Throughputmesh3D(:,:,np) = reshape(temp',[],length(temp_ue_x));
    end
    h1 = subplot(1,2,nis);
    % set(h1,'position',[0.07,0.58,0.38,0.38])
    slice(Xmesh3D,Ymesh3D,Zmesh3D,Throughputmesh3D,[],[],Hue);
    axis([min(temp_ue_x),max(temp_ue_x),min(temp_ue_y),max(temp_ue_y)])
    view(3); %grid on;
    colormap default
    hcolor=colorbar;
    set(get(hcolor,'Title'),'string','bps/Hz','FontSize',24);
    caxis([0,25])
    set(gca,'FontSize',24)
    shading interp % remove the color of frame
    xlabel("$x$(m)",'interpreter','latex','FontSize',30)
    ylabel("$y$(m)",'interpreter','latex','FontSize',30)
    zlabel("$z$(m)",'interpreter','latex','FontSize',30)
end

% head = ['Throughput(bps/Hz) distribution for IRS-aided BS with ERP:cos$\theta$'];
% head = ['cos$\theta$'];
% title(head,'FontSize',20,'interpreter','latex')

% head = ['Throughput(bps/Hz) distribution for IRS-aided BS with ERP:cos$^3\theta$'];
% head = ['cos$^3\theta$'];
% title(head,'FontSize',20,'interpreter','latex')

%% benchmarks
figure(3);
set(gcf,'Color','w')
for nbm = 1:numBenchmark
    for np = 1:numPlot
        temp = NaN*ones(length(UE_2D_raw),1);
        temp(sel_idx) = avgThroughput(:,nbm+2,np);
        Throughputmesh3D(:,:,np) = reshape(temp',[],length(temp_ue_x));
    end
    h1 = subplot(1,2,nbm);
    % set(h1,'position',[0.07,0.58,0.38,0.38])
    slice(Xmesh3D,Ymesh3D,Zmesh3D,Throughputmesh3D,[],[],Hue);
    axis([min(temp_ue_x),max(temp_ue_x),min(temp_ue_y),max(temp_ue_y)])
    view(3); %grid on;
    colormap default
    hcolor=colorbar;
    set(get(hcolor,'Title'),'string','bps/Hz','FontSize',24);
    caxis([0,25])
    set(gca,'FontSize',24)
    shading interp % remove the color of frame
    xlabel("$x$(m)",'interpreter','latex','FontSize',30)
    ylabel("$y$(m)",'interpreter','latex','FontSize',30)
    zlabel("$z$(m)",'interpreter','latex','FontSize',30)
end

% head = ['Throughput(bps/Hz) distribution for BS with fixed antenna pattern'];
% head = ['固定天线增益模式基站方案']
% title(head,'FontSize',20,'interpreter','latex')

% head = ['Throughput(bps/Hz) distribution for BS with 3D beamforming'];
% head = ['密集多天线基站三维波束成形基站方案'];
% title(head,'FontSize',20,'interpreter','latex')
