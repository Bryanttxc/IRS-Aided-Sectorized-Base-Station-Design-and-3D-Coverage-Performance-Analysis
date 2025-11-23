% created by Chen Xintong on 2022-5-3
% modified on 2022-9-25, mainly in renaming variables
% modified on 2023-07-04, mainly in renaming variables
% modified on 2023-12-17, mainly in testing code
% modified on 2025-11-22, mainly reFormat code

clear;
clc;
close all;

addpath("Geometry");
addpath("ERP");
addpath("Pathloss");

para_init;

for num_IRS_idx = 3

    fprintf("[Info] numIRS: %d x %d\n", M(num_IRS_idx), N(num_IRS_idx));

    %% 1) Generate locations of TX, IRS and UE & angles
    [UE_2D_raw, UE_2D, sel_idx, dis_irs_tx_3D, b, theta_t, theta_tx] = CalLocationAndAngle(num_IRS_idx);

    %% 2) Power gain of scattered waves
    [EA, EC] = CalPowerGainofScatteredWave(1, num_IRS_idx);

    %% 3) Ergodic Throughput
    num_UEs = size(UE_2D, 1);

    R = zeros(num_UEs,num_IRSscheme,num_trials); % instant rate
    bar_S = zeros(num_UEs,num_IRSscheme,num_trials); % mean signal power

    R_bench = zeros(num_UEs,num_benchmark,num_trials); % instantaneous rate for benchmarks
    bar_S_bench = zeros(num_UEs,num_benchmark,num_trials); % mean signal power for benchmarks

    bar_bar_S = zeros(num_UEs,num_IRSscheme+num_benchmark,num_plot); % average mean signal power
    avg_erg_throughput = zeros(num_UEs,num_IRSscheme+num_benchmark,num_plot); % average ergodic throughput

    % need 20min with 6 cpu cores
    for np = 1:num_plot
        
        tic
        fprintf("%%%%%%%%%%%%%% No. %d Layer %%%%%%%%%%%%%%\n", np);
    
        UE_z = Hue(np);
        H_u = UE_z * ones(num_UEs,1);
    
        % 3) -1] distance & angle
        % Sector represents the IRS-aided BS
        % calculate 2D, 3D distance between the sector and each UE
        Sector_2D = [0 0];
        vec_Sector2UE_2D = UE_2D-Sector_2D; %% can be negative
        vec_Sector2UE_z = H_u-H_i; %% can be negative
        dis_sector2UE_2D = sqrt(vec_Sector2UE_2D.^2);
        dis_sector2UE_3D = sqrt(sum(dis_sector2UE_2D.^2,2) + vec_Sector2UE_z.^2);
        
        % IRS-UE angle, far field, assumed to be the same for all IRS elements
        mat_Sector2UE = [vec_Sector2UE_2D vec_Sector2UE_z];
        theta_r = acosd( mat_Sector2UE*b' ./ dis_sector2UE_3D ); 
        
        % angleV, angleH of Fixed BS Pattern & 3D beamforming scheme
        angleV_atSector_local = acosd(vec_Sector2UE_z ./ dis_sector2UE_3D); % degree
        angleH_atSector_local = atan2d(vec_Sector2UE_2D(:,2), vec_Sector2UE_2D(:,1)); % degree
    
        % 3) -2] Rician K-factor
        if UE_z <= H_i  % 3GPP TR25.996 formula for ground UE
            Kfactor_LOS = 10.^((13-0.03*dis_sector2UE_3D)./10);
        else            %《3D trajectory ...》 for aerial UE
            theta_K = asin(abs(H_u-H_i)./dis_sector2UE_3D);
            Kfactor_LOS = funcRicianK(theta_K);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for gi = 1:length(G_i)
    
            % 3) -3] radiation pattern only for IRS-aided BS
            F_tx = GetRadiationValue(theta_tx, G_tx);
            F_t = GetRadiationValue(theta_t, G_i(gi));
            F_r = GetRadiationValue(theta_r, G_i(gi));
            
            EA_temp = EA(gi);
            EC_temp = EC(gi);
            parfor nt = 1:num_trials
                fprintf("No.%d trial start!!\n", nt);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 3) -4] LoS Probability & LargeScale Pathloss
                Prob_LOS = zeros(num_UEs,num_sectors);
                PL_overall = zeros(num_UEs,num_sectors);
                is_LOS = zeros(num_UEs,num_sectors);
                for ue = 1:num_UEs
                    for sec = 1:num_sectors
                        Prob_LOS(ue,sec) = LOSprobability(UE_z,dis_sector2UE_2D(ue,sec));
                        if rand(1) < Prob_LOS(ue,sec) % LOS
                            PL_overall(ue,sec) = LOSpathloss(UE_z,dis_sector2UE_2D(ue,sec),dis_sector2UE_3D(ue,sec),H_i,f_c,c);
                            is_LOS(ue,sec) = 1;
                        else % NLOS
                            PL_overall(ue,sec) = NLOSpathloss(UE_z,dis_sector2UE_2D(ue,sec),dis_sector2UE_3D(ue,sec),H_i,f_c,c);
                            is_LOS(ue,sec) = 0;
                        end
                    end
                end
                LargeScale_dB = -PL_overall;
                LargeScale = 10.^(LargeScale_dB ./ 10);
                
                idx_LoS = find(is_LOS == 1); % LOS path
                idx_NLoS = find(is_LOS == 0); % NLOS path

            %% 3) -5] -1} IRS-aided BS with cos & cos3 ERP

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % impact analysis of antenna pattern
                G_k = G_i(gi) .* F_r ./ EA_temp; % Rician K-factor gain

                Kfactor_prim = zeros(num_UEs,num_sectors); % Kfactor considering ERP
                G_rho = EA_temp.*ones(num_UEs,num_sectors); % SNR gain
                a2 = zeros(num_UEs,num_sectors);
                sigma2 = 1/2.*ones(num_UEs,num_sectors);

                G_rho(idx_LoS) = Kfactor_LOS(idx_LoS) ./ (Kfactor_LOS(idx_LoS)+1) .*G_i(gi).*F_r(idx_LoS) + 1./(Kfactor_LOS(idx_LoS)+1) .*EA_temp;
                Kfactor_prim(idx_LoS) = G_k(idx_LoS) .* Kfactor_LOS(idx_LoS);
                a2(idx_LoS) = Kfactor_prim(idx_LoS)./(2.*(Kfactor_prim(idx_LoS)+1)); % noncentral factor
                sigma2(idx_LoS) = 1./(2.*(Kfactor_prim(idx_LoS)+1)); 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generate random channel & SNR
                RcwPwr_temp = zeros(num_UEs,num_fading);
                for fad = 1:num_fading
                    fading = bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(num_UEs,M(num_IRS_idx)*N(num_IRS_idx)) ) ) + ...
                            1j* bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(num_UEs,M(num_IRS_idx)*N(num_IRS_idx)) ) );
                    RcwPwr_temp(:,fad) = P_t*G_tx*G_i(gi)*beta0*A^2.*G_rho.* LargeScale .* abs( abs(fading) * (sqrt(F_tx .* F_t) ./ dis_irs_tx_3D) ).^2;
                end
                RcwPwr = mean(RcwPwr_temp,2);
                bar_S(:,gi,nt) = RcwPwr;
                R(:,gi,nt) = log2(1 + RcwPwr./NoisePower);

            %% 3) -5] -2} Fixed BS Pattern & 3D beamforming
                is_FixedPattern = 0;
                is_3Dbeamforming = 0;

                % Sector antennas gain
                if gi+2 == 3 % Fixed Pattern
                    is_FixedPattern = 1;
                    GsectordB = ArrayPowerPattern(angleV_atSector_local,angleH_atSector_local,angleTiltV,angleTiltH,M_b(num_IRS_idx),N_b(num_IRS_idx),dV,dH,wavelength);
                elseif gi+2 == 4 % 3D beamforming
                    is_3Dbeamforming = 1;
                    GsectordB = ElementPowerPatternOverall(angleV_atSector_local,angleH_atSector_local); % the element power pattern
                end

                % calculate the parameters accounting for the impact of radiation pattern
                G_k = 10.^(GsectordB./10) ./ EC_temp; % Ricial K-factor gain

                G_rho       = zeros(num_UEs,num_sectors);
                Kfactor_prim = zeros(num_UEs,num_sectors);
                a2          = zeros(num_UEs,num_sectors);
                sigma2      = 1/2.*ones(num_UEs,num_sectors);

                if is_FixedPattern
                    G_rho(idx_LoS) = ( Kfactor_LOS(idx_LoS)./(Kfactor_LOS(idx_LoS)+1) .*(10.^(GsectordB(idx_LoS)./10)) + 1./(Kfactor_LOS(idx_LoS)+1) .*EC_temp ) ./ (M_b(num_IRS_idx)*N_b(num_IRS_idx)); % SNR gain
                    G_rho(idx_NLoS) = EC_temp / (M_b(num_IRS_idx)*N_b(num_IRS_idx)); % SNR gain
                elseif is_3Dbeamforming
                    G_rho(idx_LoS) = Kfactor_LOS(idx_LoS)./(Kfactor_LOS(idx_LoS)+1) .*(10.^(GsectordB(idx_LoS)./10)) + 1./(Kfactor_LOS(idx_LoS)+1) .*EC_temp; % SNR gain
                    G_rho(idx_NLoS) = EC_temp; % SNR gain
                end
                Kfactor_prim(idx_LoS) = G_k(idx_LoS) .* Kfactor_LOS(idx_LoS);
                a2(idx_LoS) = Kfactor_prim(idx_LoS)./(2.*(Kfactor_prim(idx_LoS)+1)); % noncentral factor
                sigma2(idx_LoS) = 1./(2.*(Kfactor_prim(idx_LoS)+1));

                % calculate each UE SNR
                if is_FixedPattern
                    fading = ( sqrt(a2)+sqrt(sigma2).*randn(num_UEs,num_fading) ) + 1j* ( sqrt(a2)+sqrt(sigma2).*randn(num_UEs,num_fading) );
                    RcwPwr_temp = P_t .* G_rho .* LargeScale .* abs(fading).^2; % numUEs x num_fading
                elseif is_3Dbeamforming
                   for fad = 1:num_fading
                       ch_3D_temp = sqrt(LargeScale) .* sqrt(G_rho) .* ( bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(num_UEs,M_b(num_IRS_idx)*N_b(num_IRS_idx)) )) +...
                                             1j* bsxfun( @plus,sqrt(a2),bsxfun( @times,sqrt(sigma2),randn(num_UEs,M_b(num_IRS_idx)*N_b(num_IRS_idx)) )) );
                       RcwPwr_temp(:,fad) = P_t .* sum(abs(ch_3D_temp).^2, 2); % MRT beamforming, norm(Hi).^2
                   end
                end
                RcvPwr = mean(RcwPwr_temp,2);
                bar_S_bench(:,gi,nt) = RcvPwr;
                R_bench(:,gi,nt) = log2(1 + RcvPwr./NoisePower);
            end
        end

        % 3) -6] ergodic throughput
        avg_erg_throughput(:,1:num_IRSscheme,np) = mean(R,3);
        avg_erg_throughput(:,num_IRSscheme+1:num_IRSscheme+2,np) = mean(R_bench,3);
        bar_bar_S(:,1:num_IRSscheme,np) = mean(bar_S,3);
        bar_bar_S(:,num_IRSscheme+1:num_IRSscheme+2,np) = mean(bar_S_bench,3);
        toc
    end

end

%% 4) save data
date = char(datetime('today', 'InputFormat','yyyy-MM-dd'));
save(['../data/', date, '_avgThroughputinSingleCell'],'avg_erg_throughput');
save(['../data/', date, '_barbarSinSingleCell'],'bar_bar_S');

%% 5) performance evaluation
% Average Ergodic Throughput (AET) bps/Hz
AETRate(1,:) = sum(avg_erg_throughput(:,:,1)) / num_UEs; % Ground UE
AETRate(2,:) = sum(avg_erg_throughput(:,:,2)) / num_UEs; % Aerial UE

% Jains Fairness
Fairness(1,:) = sum(avg_erg_throughput(:,:,1)).^2 ./ ( length(avg_erg_throughput) .* sum(avg_erg_throughput(:,:,1).^2) ); % Ground UE
Fairness(2,:) = sum(avg_erg_throughput(:,:,2)).^2 ./ ( length(avg_erg_throughput) .* sum(avg_erg_throughput(:,:,2).^2) ); % Aerial UE

% Average Mean Signal Power (dB)
AvgMeanSignalPower(1,:) = 10.*log10( mean(bar_bar_S(:,:,1),1) ); % Ground UE
AvgMeanSignalPower(2,:) = 10.*log10( mean(bar_bar_S(:,:,2),1) ); % Aerial UE

AETRate
Fairness
AvgMeanSignalPower

% Ergodic Throughput distribution plot
temp_ue_x = 0:3:2*r; % x
temp_ue_y = -sqrt(3)/2*r:10:sqrt(3)/2*r; % y

%%% 5) -1] Fixed Pattern
figure(2);
set(gcf,'Color','w')
[Xmesh3D,Ymesh3D,Zmesh3D] = meshgrid(temp_ue_x,temp_ue_y,Hue);

for np = 1:num_plot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avg_erg_throughput(:,3,np);
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

%%% 5) -2] 3D beamforming
for np = 1:num_plot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avg_erg_throughput(:,4,np);
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

%%% 5) -3] IRS-aided BS with cos theta
figure(3);
set(gcf,'Color','w')
for np = 1:num_plot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avg_erg_throughput(:,1,np);
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

%%% 5) -4] IRS-aided BS with cos3 theta
for np = 1:num_plot
    temp = NaN*ones(length(UE_2D_raw),1);
    temp(sel_idx) = avg_erg_throughput(:,2,np);
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
