% created by Chen Xintong on 2022-5-3
% modified on 2022-9-25, mainly in renaming variables
% modified on 2023-07-04, mainly in renaming variables
% modified on 2023-12-17, mainly in testing code
clc;
clear;
close all

para_init;

for nnn = 3
nnn

%% Generate location of TX,IRS and UE & angle
[UE_2D_raw, UE_2D, sel_idx, dis_irs_tx_3D, b, theta_t, theta_tx] = calLocationAndAngle(nnn);

%% Power gain of scattered waves
[EA, EC] = calPowerGainofScatteredWave(1, nnn);

%% Ergodic Throughput
numUEs = size(UE_2D,1); % number of UEs

R = zeros(numUEs,numIRSscheme,numTrials); % instantaneous rate
barS = zeros(numUEs,numIRSscheme,numTrials); % mean signal power

R_bench = zeros(numUEs,numBenchmark,numTrials); % instantaneous rate for benchmarks
barS_bench = zeros(numUEs,numBenchmark,numTrials); % mean signal power for benchmarks

barbarS = zeros(numUEs,numIRSscheme+numBenchmark,numPlot); % average mean signal power
avgThroughput = zeros(numUEs,numIRSscheme+numBenchmark,numPlot); % average ergodic throughput

% need 20min with 6 cores
for np = 1:numPlot
    
    tic
    fprintf("%%%%%%%%%%%%%% No. %d Layer %%%%%%%%%%%%%%\n",np)

    UE_z = Hue(np);
    H_u = UE_z * ones(numUEs,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sector represents the IRS-aided BS
    % calculate 2D, 3D distance between the sector and each UE
    Sector_2D = [0 0];
    vec_SectortoUE_2D = UE_2D-Sector_2D; %% can be negative
    vec_SectortoUE_z = H_u-H_i; %% can be negative
    dis_sector_ue_2D = sqrt(vec_SectortoUE_2D.^2);
    dis_sector_ue_3D = sqrt(sum(dis_sector_ue_2D.^2,2)+vec_SectortoUE_z.^2);
    
    % IRS-UE angle, far field, assumed to be the same for all IRS elements
    mat_SectorToUE = [vec_SectortoUE_2D vec_SectortoUE_z];
    theta_r = acosd( mat_SectorToUE*b' ./ dis_sector_ue_3D ); 
    
    % angleV, angleH of Fixed BS Pattern & 3D beamforming scheme
    angleVatSectorLocal = acosd(vec_SectortoUE_z./dis_sector_ue_3D); % degree
    angleHatSectorLocal = atan2d(vec_SectortoUE_2D(:,2),vec_SectortoUE_2D(:,1)); % degree

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rician K-factor
    if UE_z <= H_i  % 3GPP TR25.996 formula for ground UE
        Kfactor_LOS = 10.^((13-0.03*dis_sector_ue_3D)./10);
    else            %《3D trajectory ...》 for aerial UE
        theta_K = asin(abs(H_u-H_i)./dis_sector_ue_3D);
        Kfactor_LOS = funcRicianK(theta_K);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for gi = 1:length(G_i)

        % radiation pattern
        F_tx = getRadiationValue(theta_tx, G_tx);
        F_t = getRadiationValue(theta_t, G_i(gi));
        F_r = getRadiationValue(theta_r, G_i(gi));
        
        EA_temp = EA(gi);
        EC_temp = EC(gi);
        parfor nt = 1:numTrials
            nt
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LoS Probability & LargeScale Pathloss
            PrLOS = zeros(numUEs,numSectors); % LOS Probability
            PL_overall = zeros(numUEs,numSectors); % Pathloss
            isLOS = zeros(numUEs,numSectors);
            for p = 1:numUEs
                for q = 1:numSectors
                    PrLOS(p,q) = LOSprobability(UE_z,dis_sector_ue_2D(p,q));
                    if rand(1) < PrLOS(p,q) %% LOS
                        PL_overall(p,q) = LOSpathloss(UE_z,dis_sector_ue_2D(p,q),dis_sector_ue_3D(p,q),H_i,f_c,c);
                        isLOS(p,q) = 1;
                    else %% NLOS
                        PL_overall(p,q) = NLOSpathloss(UE_z,dis_sector_ue_2D(p,q),dis_sector_ue_3D(p,q),H_i,f_c,c);
                        isLOS(p,q) = 0;
                    end
                end
            end
            LargeScaledB = -PL_overall;
            LargeScale = 10.^(LargeScaledB ./ 10);
            
            idxLoS = find(isLOS == 1); % LOS path
            idxNLoS = find(isLOS == 0); % NLOS path
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % impact analysis of antenna pattern
            G_k = G_i(gi) .* F_r ./ EA_temp; % Rician K-factor gain

            KfactorPrim = zeros(numUEs,numSectors); % Kfactor considering radiation pattern
            G_rho = EA_temp.*ones(numUEs,numSectors); % SNR gain
            a2 = zeros(numUEs,numSectors);
            sigma2 = 1/2.*ones(numUEs,numSectors);

            G_rho(idxLoS) = Kfactor_LOS(idxLoS) ./ (Kfactor_LOS(idxLoS)+1) .*G_i(gi).*F_r(idxLoS) + 1./(Kfactor_LOS(idxLoS)+1) .*EA_temp;
            KfactorPrim(idxLoS) = G_k(idxLoS) .* Kfactor_LOS(idxLoS);
            a2(idxLoS) = KfactorPrim(idxLoS)./(2.*(KfactorPrim(idxLoS)+1)); % noncentral factor
            sigma2(idxLoS) = 1./(2.*(KfactorPrim(idxLoS)+1)); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate random channel & SNR
            RcwPwr_temp = zeros(numUEs,numfading);
            for ran = 1:numfading
                fading = bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(numUEs,M(nnn)*N(nnn)) ) ) + ...
                        1j* bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(numUEs,M(nnn)*N(nnn)) ) );
                RcwPwr_temp(:,ran) = P_t*G_tx*G_i(gi)*beta0*A^2.*G_rho.* LargeScale .* abs( abs(fading) * (sqrt(F_tx .* F_t) ./ dis_irs_tx_3D) ).^2;
            end
            RcwPwr = mean(RcwPwr_temp,2);
            barS(:,gi,nt) = RcwPwr;
            R(:,gi,nt) = log2(1 + RcwPwr./NoisePower);

       %% Fixed BS Pattern & 3D beamforming
            isFixedPattern = 0;
            is3Dbeamforming = 0;
            
            % Sector antennas' gain
            if gi+2 == 3 %% Fixed Pattern
                isFixedPattern = 1;
                GsectordB = ArrayPowerPattern(angleVatSectorLocal,angleHatSectorLocal,angleTiltV,angleTiltH,M_b(nnn),N_b(nnn),dV,dH,wavelength);
            elseif gi+2 == 4 %% 3D beamforming
                is3Dbeamforming = 1;
                GsectordB = ElementPowerPatternOverall(angleVatSectorLocal,angleHatSectorLocal); % the element power pattern
            end

            % calculate the parameters accounting for the impact of radiation pattern
            G_k = 10.^(GsectordB./10) ./ EC_temp; % Ricial K-factor gain
            
            G_rho       = zeros(numUEs,numSectors);
            KfactorPrim = zeros(numUEs,numSectors);
            a2          = zeros(numUEs,numSectors);
            sigma2      = 1/2.*ones(numUEs,numSectors);
            
            if isFixedPattern
                G_rho(idxLoS) = ( Kfactor_LOS(idxLoS)./(Kfactor_LOS(idxLoS)+1) .*(10.^(GsectordB(idxLoS)./10)) + 1./(Kfactor_LOS(idxLoS)+1) .*EC_temp ) ./ (M_b(nnn)*N_b(nnn)); % SNR gain
                G_rho(idxNLoS) = EC_temp / (M_b(nnn)*N_b(nnn)); % SNR gain
            elseif is3Dbeamforming
                G_rho(idxLoS) = Kfactor_LOS(idxLoS)./(Kfactor_LOS(idxLoS)+1) .*(10.^(GsectordB(idxLoS)./10)) + 1./(Kfactor_LOS(idxLoS)+1) .*EC_temp; % SNR gain
                G_rho(idxNLoS) = EC_temp; % SNR gain
            end
            KfactorPrim(idxLoS) = G_k(idxLoS) .* Kfactor_LOS(idxLoS);
            a2(idxLoS) = KfactorPrim(idxLoS)./(2.*(KfactorPrim(idxLoS)+1)); % noncentral factor
            sigma2(idxLoS) = 1./(2.*(KfactorPrim(idxLoS)+1));

            % calculate each UE's SNR
            if isFixedPattern
                fading = ( sqrt(a2)+sqrt(sigma2).*randn(numUEs,numfading) ) + 1j* ( sqrt(a2)+sqrt(sigma2).*randn(numUEs,numfading) );
                RcwPwr_temp = P_t .* G_rho .* LargeScale .* abs(fading).^2; % numUEs x numFading
            elseif is3Dbeamforming
               for ran = 1:numfading
                   ch_3D_temp = sqrt(LargeScale) .* sqrt(G_rho) .* ( bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(numUEs,M_b(nnn)*N_b(nnn)) )) +...
                                         1j* bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(numUEs,M_b(nnn)*N_b(nnn)) )) );
                   RcwPwr_temp(:,ran) = P_t .* sum(abs(ch_3D_temp).^2, 2); % MRT beamforming, norm(Hi).^2
               end
            end
            RcvPwr = mean(RcwPwr_temp,2);
            barS_bench(:,gi,nt) = RcvPwr;
            R_bench(:,gi,nt) = log2(1 + RcvPwr./NoisePower);
        end
    end
    
    avgThroughput(:,1:numIRSscheme,np) = mean(R,3);
    avgThroughput(:,numIRSscheme+1:numIRSscheme+2,np) = mean(R_bench,3);
    barbarS(:,1:numIRSscheme,np) = mean(barS,3);
    barbarS(:,numIRSscheme+1:numIRSscheme+2,np) = mean(barS_bench,3);
    toc
end

end

%% save
save('../data/20231217avgThroughputinSingleCell','avgThroughput');
save('../data/20231217barbarSinSingleCell','barbarS');

%% performance evaluation
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

figure(2);
set(gcf,'Color','w')
[Xmesh3D,Ymesh3D,Zmesh3D] = meshgrid(temp_ue_x,temp_ue_y,Hue);

for np = 1:numPlot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avgThroughput(:,3,np);
    Throughputmesh3D(:,:,np) = reshape(temp',[],length(temp_ue_x));
end
h1 = subplot(1,2,1);
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
% head = ['Throughput(bps/Hz) distribution for BS with fixed antenna pattern'];
% head = ['固定天线增益模式基站方案']
% title(head,'FontSize',20,'interpreter','latex')

for np = 1:numPlot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avgThroughput(:,4,np);
    Throughputmesh3D(:,:,np) = reshape(temp',[],length(temp_ue_x));
end
h2 = subplot(1,2,2);
% set(h2,'position',[0.55,0.58,0.38,0.38])
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
% head = ['Throughput(bps/Hz) distribution for BS with 3D beamforming'];
% head = ['密集多天线基站三维波束成形基站方案'];
% title(head,'FontSize',20,'interpreter','latex')

figure(3);
set(gcf,'Color','w')
for np = 1:numPlot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avgThroughput(:,1,np);
    Throughputmesh3D(:,:,np) = reshape(temp',[],length(temp_ue_x));
end
h3 = subplot(1,2,1);
% set(h3,'position',[0.07,0.0875,0.38,0.38])
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
% head = ['Throughput(bps/Hz) distribution for IRS-mounted BS with ERP:cos$\theta$'];
% head = ['cos$\theta$'];
% title(head,'FontSize',20,'interpreter','latex')

for np = 1:numPlot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avgThroughput(:,2,np);
    Throughputmesh3D(:,:,np) = reshape(temp',[],length(temp_ue_x));
end
h4 = subplot(1,2,2);
% set(h4,'position',[0.55,0.0875,0.38,0.38])
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
% head = ['Throughput(bps/Hz) distribution for IRS-mounted BS with ERP:cos$^3\theta$'];
% head = ['cos$^3\theta$'];
% title(head,'FontSize',20,'interpreter','latex')

