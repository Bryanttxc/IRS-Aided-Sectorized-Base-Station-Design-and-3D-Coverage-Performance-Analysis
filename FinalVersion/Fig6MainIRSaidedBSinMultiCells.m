%created by Chen Xintong on 2023-12-17
% modified on 2025-11-22, mainly reFormat code

clear;
clc;
close all;

addpath("Geometry");
addpath("ERP");
addpath("Pathloss");

para_init;

%% 1) Generate sector locations in multiple cells
[sector_pos, label_pos, num_sectors] = GenerateSectorLoc();

%% 2) Generate locations of TX,IRS and UE & angles
num_IRS_idx = 3;
[UE_2D_raw, UE_2D, sel_idx, dis_irs_tx_3D, b, theta_t, theta_tx] = CalLocationAndAngle(num_IRS_idx);

% boresight of each sector
% b1_i : sector 1 4 7 10 13 16 19
% b2_i : sector 2 5 8 11 14 17 20
% b3_i : sector 3 6 9 12 15 18 21

b1_i = [cosd(0)*cosd(alpha_y)   sind(0)*cosd(alpha_y)   sind(alpha_y)];
b2_i = [cosd(120)*cosd(alpha_y) sind(120)*cosd(alpha_y) sind(alpha_y)];
b3_i = [cosd(240)*cosd(alpha_y) sind(240)*cosd(alpha_y) sind(alpha_y)];

%% 3) Power gain of scattered waves
[EA, EC] = CalPowerGainofScatteredWave(1, num_IRS_idx);

%% 4) Ergodic Throughput
num_row = 2; % Row1: the channel between the UE in cell 1 and all sectors, Row2: the channel between each sector and its own served UE 
num_UEs = num_sectors; % assume each Sector have one UE
R_temp = zeros(num_trials,num_plot); % average Throughput
R_temp_bench = zeros(num_trials,num_plot); % average Throughput
R = zeros(num_trials,num_plot,num_IRSscheme+num_benchmark); % Four scheme
avgThroughput = zeros(size(UE_2D,1),num_IRSscheme+num_benchmark,num_plot);

barS_irs_temp = zeros(num_trials,num_plot);
barI_irs_temp = zeros(num_trials,num_plot);
barS_bench_temp = zeros(num_trials,num_plot);
barI_bench_temp = zeros(num_trials,num_plot);

barS = zeros(num_trials,num_plot,num_IRSscheme+num_benchmark); % Four scheme
barI = zeros(num_trials,num_plot,num_IRSscheme+num_benchmark);

barbarS = zeros(size(UE_2D,1),num_IRSscheme+num_benchmark,num_plot);
barbarI = zeros(size(UE_2D,1),num_IRSscheme+num_benchmark,num_plot);

mean_sigPow_lb_temp = zeros(num_trials,num_plot); % lower bound
mean_itfPow_ub_temp = zeros(num_trials,num_plot); % upper bound
mean_sigPow_lb = zeros(num_trials,num_plot,num_IRSscheme+num_benchmark); % Four scheme
mean_itfPow_ub = zeros(num_trials,num_plot,num_IRSscheme+num_benchmark);
avg_mean_SigPow_lb = zeros(size(UE_2D,1),num_IRSscheme+num_benchmark,num_plot);
avg_mean_ItrPow_ub = zeros(size(UE_2D,1),num_IRSscheme+num_benchmark,num_plot);

% need very long time
for ue = 1:size(UE_2D,1)

    tic
    fprintf("%%%%%%%%%%%%%% No. %d UE %%%%%%%%%%%%%%\n", ue);
    
    % 4) -1] Generate UE random locations of other Sectors
    UE_max_radius = r * sqrt(3)/2; % inscribed circle
    theta = rand(num_UEs-1,num_trials)* 2*pi; % uniform distribution within the circle
    r_ue = rand(num_UEs-1,num_trials).* UE_max_radius;
    x_ue = r_ue.*cos(theta) + label_pos(2:end,1);
    y_ue = r_ue.*sin(theta) + label_pos(2:end,2);

    otherSector_UE_point(:,:,1) = x_ue;
    otherSector_UE_point(:,:,2) = y_ue;
    
    sel_UE_point(:,:,1) = [repmat(UE_2D(ue,1),1,num_trials);otherSector_UE_point(:,:,1)];
    sel_UE_point(:,:,2) = [repmat(UE_2D(ue,2),1,num_trials);otherSector_UE_point(:,:,2)];
    
    for gi = 1:length(G_i)
        G_i_temp = G_i(gi);
        EA_temp = EA(gi);
        EC_temp = EC(gi);
        
        % 4) -2] radiation pattern only for IRS-aided BS
        F_tx = GetRadiationValue(theta_tx, G_tx);
        F_t  = GetRadiationValue(theta_t, G_i_temp);
        [MinFcombine,Minidx] = min(F_tx.*F_t);
        [MaxFcombine,Maxidx] = max(F_tx.*F_t);

        parfor nt = 1:num_trials

            % 4) -3] distances & angles
            % 2D,3D distance and angleH, angleV from each UE to each sector
            vec_sector2UE_x = bsxfun(@minus, sel_UE_point(:,nt,1), sector_pos(:,1)'); %% can be negative
            vec_sector2UE_y = bsxfun(@minus, sel_UE_point(:,nt,2), sector_pos(:,2)'); %% can be negative

            % delete unnecessary data
            sel_vec_sector2UE_x = [vec_sector2UE_x(1,:);diag(vec_sector2UE_x)'];
            sel_vec_sector2UE_y = [vec_sector2UE_y(1,:);diag(vec_sector2UE_y)'];

            dist_sector2UE_2D = sqrt(sel_vec_sector2UE_x.^2+sel_vec_sector2UE_y.^2);

            % angleH of Fixed BS Pattern & 3D beamforming
            angleH_atSector_global = atan2d(sel_vec_sector2UE_y,sel_vec_sector2UE_x); % degree
            angleH_atSector_local = angleH_atSector_global; %% relative to the boresight direction of each sector, 0,120,-120
            angleH_atSector_local(:,2:3:end) = angleH_atSector_global(:,2:3:end)-120;
            angleH_atSector_local(:,3:3:end) = angleH_atSector_global(:,3:3:end)-(-120);

            % convert to -180~180
            idx_negt = find(angleH_atSector_local<-180);
            angleH_atSector_local(idx_negt) = angleH_atSector_local(idx_negt)+360;
            idx_greater180 = find(angleH_atSector_local>180);
            angleH_atSector_local(idx_greater180) = angleH_atSector_local(idx_greater180)-360;
            index = find(angleH_atSector_local < -90 | angleH_atSector_local > 90);

            for np = 1:num_plot

                UE_z = Hue(np)*ones(num_UEs,1);
                vec_sector2UE_z = bsxfun(@minus, UE_z, H_i*ones(num_sectors,1)'); %% can be negative
                sel_vec_sector2UE_z = [vec_sector2UE_z(1,:);diag(vec_sector2UE_z)'];
                dist_sector2UE_3D = sqrt(dist_sector2UE_2D.^2+sel_vec_sector2UE_z.^2);

                %%% The angleV of FixedPattern & 3D beamforming BS %%%
                angleV_atSector_local=acosd(sel_vec_sector2UE_z./dist_sector2UE_3D);

                %%%%%%%%%%%%%%%%%%%% The angleV of cos & cos3 IRS %%%%%%%%%%%%%%%%%%%%
                sector1ToUE = [reshape(sel_vec_sector2UE_x(:,1:3:end),[],1) reshape(sel_vec_sector2UE_y(:,1:3:end),[],1) reshape(sel_vec_sector2UE_z(:,1:3:end),[],1)];
                sector2ToUE = [reshape(sel_vec_sector2UE_x(:,2:3:end),[],1) reshape(sel_vec_sector2UE_y(:,2:3:end),[],1) reshape(sel_vec_sector2UE_z(:,2:3:end),[],1)]; 
                sector3ToUE = [reshape(sel_vec_sector2UE_x(:,3:3:end),[],1) reshape(sel_vec_sector2UE_y(:,3:3:end),[],1) reshape(sel_vec_sector2UE_z(:,3:3:end),[],1)];

                theta_r = zeros(num_row,num_sectors);
                theta_r(:,1:3:end) = reshape( acosd( sector1ToUE*b1_i' ./ reshape(dist_sector2UE_3D(:,1:3:end),[],1) )' ,num_row,[]);
                theta_r(:,2:3:end) = reshape( acosd( sector2ToUE*b2_i' ./ reshape(dist_sector2UE_3D(:,2:3:end),[],1) )' ,num_row,[]);
                theta_r(:,3:3:end) = reshape( acosd( sector3ToUE*b3_i' ./ reshape(dist_sector2UE_3D(:,3:3:end),[],1) )' ,num_row,[]);

                % 4) -4] LOS Probability and Large-Scale Pathloss
                Prob_LOS = zeros(num_row,num_sectors);
                PL_overall = zeros(num_row,num_sectors);
                isLOS = zeros(num_row,num_sectors);
                for p = 1:num_row
                    for q = 1:num_sectors
                        Prob_LOS(p,q) = LOSprobability(UE_z(1),dist_sector2UE_2D(p,q));
                    end
                end
                Prob_LOS(index) = 0; % modify
                
                for p = 1:num_row
                    for q = 1:num_sectors
                        if rand(1) < Prob_LOS(p,q) %% LOS
                            PL_overall(p,q) = LOSpathloss(UE_z(1),dist_sector2UE_2D(p,q),dist_sector2UE_3D(p,q),H_i,f_c,c);
                            isLOS(p,q) = 1;
                        else %% NLOS
                            PL_overall(p,q) = NLOSpathloss(UE_z(1),dist_sector2UE_2D(p,q),dist_sector2UE_3D(p,q),H_i,f_c,c);
                            isLOS(p,q) = 0;
                        end
                    end
                end
                LargeScaledB = -PL_overall;
                LargeScale = 10 .^(LargeScaledB ./ 10);

                % 4) -5] Rician Kfactor
                if UE_z(1) <= H_i
                    Kfactor_LOS = 10.^((13-0.03*dist_sector2UE_3D)./10);
                else
                    theta_K = asin(abs(UE_z(1)'-H_i)./dist_sector2UE_3D);
                    Kfactor_LOS = funcRicianK(theta_K);
                end
                Kfactor_LOS(index) = 0;
                Kfactor_NLOS = zeros(num_row,num_sectors); 

            %% 4) -6] -1} IRS-aided BS

                % IRS cos & cos3 ERP only for IRS-aided BS
                F_r = GetRadiationValue(theta_r, G_i_temp);

                % Rician K-factor gain
                Gk = G_i_temp .* F_r ./ EA_temp; 

                % determine parameters of fading
                KfactorPrim = zeros(num_row,num_sectors); % Kfactor considering ERP
                Grho = zeros(num_row,num_sectors); % SNR gain
                a2 = zeros(num_row,num_sectors);
                sigma2 = zeros(num_row,num_sectors);
                for p = 1
                    for q = 1:num_sectors
                        if isLOS(p,q) %% LOS
                            Grho(p,q) = Kfactor_LOS(p,q)./(Kfactor_LOS(p,q)+1) .*G_i(gi).*F_r(p,q) + 1./(Kfactor_LOS(p,q)+1) .*EA_temp;
                            KfactorPrim(p,q) = Gk(p,q) .* Kfactor_LOS(p,q);
                            a2(p,q) = KfactorPrim(p,q)./( 2.*(KfactorPrim(p,q)+1) ); % noncentral factor
                            sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                        else %%NLOS
                            Grho(p,q) = EA_temp;
                            KfactorPrim(p,q) = Kfactor_NLOS(p,q);
                            a2(p,q) = KfactorPrim(p,q)./( 2.*(KfactorPrim(p,q)+1) ); % noncentral factor
                            sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                        end
                    end
                end

                % Generate random channel & SINR
                sigPower_temp = zeros(num_fading,1);
                intfPower_temp = zeros(num_fading,1);
                for ran = 1:num_fading
                    for p = 1
                        fading = ( sqrt(a2(p,p))+sqrt(sigma2(p,p)).*randn(M(num_IRS_idx)*N(num_IRS_idx),1) )+...
                            1j*( sqrt(a2(p,p))+sqrt(sigma2(p,p)).*randn(M(num_IRS_idx)*N(num_IRS_idx),1) );% small scale fading
                        sigPower_temp(ran,p) =  P_t * G_tx*G_i(gi)*beta0*A^2.*Grho(p,p) .* LargeScale(p,p).*...
                                      abs( sum( abs(fading).* sqrt(F_tx.*F_t) ./ dis_irs_tx_3D ) ).^2;
                        for q = [1:p-1 p+1:num_UEs] % for other UEs
                            randPhase = exp(1j*2*pi*rand(M(num_IRS_idx)*N(num_IRS_idx),1));
                            fading = ( sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M(num_IRS_idx)*N(num_IRS_idx),1) )+...
                                1j*( sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M(num_IRS_idx)*N(num_IRS_idx),1) );% small scale fading
                            intfPower_temp(ran,p) = intfPower_temp(ran,p) + P_t * G_tx*G_i(gi)*beta0.*Grho(p,q).* LargeScale(p,q).*...
                                  abs( sum( abs(fading).* sqrt(F_tx.*F_t).*A.*randPhase.*exp(-1j*2*pi/wavelength.*dis_irs_tx_3D).*exp(-1j*2*pi/wavelength.*dist_sector2UE_3D(p,q))./dis_irs_tx_3D ) ).^2;
                        end
                    end
                end
                SINR = mean( sigPower_temp ./ (intfPower_temp + NoisePower) );
                barS_irs_temp(nt,np) = mean(sigPower_temp);
                barI_irs_temp(nt,np) = mean(intfPower_temp);
                mean_sigPow_lb_temp(nt,np) = ( pi / (4*(KfactorPrim(1,1)+1))*L(-KfactorPrim(1,1))^2*(M(num_IRS_idx)*N(num_IRS_idx))^2 + (1 - pi/(4*(KfactorPrim(1,1)+1)))*L(-KfactorPrim(1,1))^2*M(num_IRS_idx)*N(num_IRS_idx) ).*...
                                        P_t* (A^2*G_tx*G_i(gi)*Grho(1,1)*MinFcombine*beta0/dis_irs_tx_3D(Minidx)*LargeScale(1,1));
                mean_itfPow_ub_temp(nt,np) = sum( M(num_IRS_idx)*N(num_IRS_idx)* P_t* (A^2*G_tx*G_i(gi)*MaxFcombine*beta0/dis_irs_tx_3D(Maxidx).*Grho(1,2:end).*LargeScale(1,2:end)) );
                R_temp(nt,np) = log2(1+SINR)';

           %% 4) -6] -2} Fixed BS Pattern & 3D beamforming

                is_FixedPattern = 0;
                is_3Dbeamforming = 0;

                % Sector antennas' Gain
                if gi+2 == 3 % FixedPattern
                    is_FixedPattern = 1;
                    Gsector_dB = zeros(num_row,num_sectors);
                    for p = 1:num_row
                        for q = 1:num_sectors
                            Gsector_dB(p,q) = ArrayPowerPattern(angleV_atSector_local(p,q),angleH_atSector_local(p,q),angleTiltV,angleTiltH,M_b(num_IRS_idx),N_b(num_IRS_idx),dV,dH,wavelength);
                        end
                    end
                elseif gi+2 == 4 % 3D beamforming
                    is_3Dbeamforming = 1;
                    Gsector_dB = zeros(num_row,num_sectors);
                    for p = 1:num_row
                        for q = 1:num_sectors
                            Gelement_dB = ElementPowerPatternOverall(angleV_atSector_local(p,q),angleH_atSector_local(p,q)); % the element power gain
                            Gelement = 10^(Gelement_dB/10);
                            Felement = sqrt(Gelement); % the element field pattern
                            rMatrix = zeros(M_b(num_IRS_idx)*N_b(num_IRS_idx),3);
                            idx = -N_b(num_IRS_idx)/2+1:N_b(num_IRS_idx)/2;
                            cnt = 1;
                            for v = 1:N_b(num_IRS_idx)
                                 rMatrix(1+(cnt-1)*M_b(num_IRS_idx):cnt*M_b(num_IRS_idx),3) = [-M_b(num_IRS_idx)/2:(M_b(num_IRS_idx)/2-1)]'*dV;
                                 rMatrix(1+(cnt-1)*M_b(num_IRS_idx):cnt*M_b(num_IRS_idx),2) = (idx(v)-1)*dH;
                                 cnt = cnt + 1;
                            end
                            k_vector = (2*pi/wavelength)*[sind(angleV_atSector_local(p,q))*cosd(angleH_atSector_local(p,q)),sind(angleV_atSector_local(p,q))*sind(angleH_atSector_local(p,q)),cosd(angleV_atSector_local(p,q))];%wave vector
                            steering_vector = exp(-1i.*(rMatrix*k_vector')); % Num_Ant x 1
                            Gsector_dB(p,q) = 10*log10( abs(Felement*steering_vector(centeridx(num_IRS_idx)))^2 );
                        end
                    end
                end

                Gk = 10.^(Gsector_dB./10) ./ EC_temp; % Ricial K-factor gain
                Gk(index) = 0;

                % determine parameters of each UE
                Grho        = zeros(num_row,num_sectors);
                KfactorPrim = zeros(num_row,num_sectors);
                a2          = zeros(num_row,num_sectors);
                sigma2      = zeros(num_row,num_sectors);
                fading      = zeros(num_row,num_sectors);
                ch_3D_temp  = zeros(M_b(num_IRS_idx)*N_b(num_IRS_idx),num_sectors,num_row,num_fading);                
                for p = 1:num_row
                    for q = 1:num_sectors
                        if isLOS(p,q)
                            if is_FixedPattern
                                Grho(p,q) = ( Kfactor_LOS(p,q)./(Kfactor_LOS(p,q)+1) .*(10.^(Gsector_dB(p,q)./10)) + 1./(Kfactor_LOS(p,q)+1) .*EC_temp ) ./ (M_b(num_IRS_idx)*N_b(num_IRS_idx)); % SNR gain
                            elseif is_3Dbeamforming
                                Grho(p,q) = Kfactor_LOS(p,q)./(Kfactor_LOS(p,q)+1) .*(10.^(Gsector_dB(p,q)./10)) + 1./(Kfactor_LOS(p,q)+1) .*EC_temp; % SNR gain
                            end
                            KfactorPrim(p,q) = Gk(p,q) .* Kfactor_LOS(p,q);
                            a2(p,q) = KfactorPrim(p,q)./(2.*(KfactorPrim(p,q)+1)); % noncentral factor
                            sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                        else
                            if is_FixedPattern
                                Grho(p,q) = EC_temp / (M_b(num_IRS_idx)*N_b(num_IRS_idx)); % SNR gain
                            elseif is_3Dbeamforming
                                Grho(p,q) = EC_temp;
                            end
                            KfactorPrim(p,q) = Gk(p,q) .* Kfactor_NLOS(p,q);                            
                            a2(p,q) = KfactorPrim(p,q)./(2.*(KfactorPrim(p,q)+1)); % noncentral factor
                            sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                        end

                        if is_3Dbeamforming
                            for ran = 1:num_fading
                                 ch_3D_temp(:,q,p,ran) = sqrt(Grho(p,q)) .* sqrt(LargeScale(p,q)) .* ...
                                     ( (sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M_b(num_IRS_idx)*N_b(num_IRS_idx),1)) + 1j* (sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M_b(num_IRS_idx)*N_b(num_IRS_idx),1)) );
                            end
                        end

                    end
                end

                % calculate each UE's SINR
                power_temp = zeros(num_row,num_sectors);
                sigPower_temp = zeros(num_fading,1);
                intfPower_temp = zeros(num_fading,1);
                if is_FixedPattern
                    for ran = 1:num_fading
                        fading = ( sqrt(a2)+sqrt(sigma2).*randn(num_row,num_sectors) ) +...
                                    1j* ( sqrt(a2)+sqrt(sigma2).*randn(num_row,num_sectors) );
                        power_temp = P_t .* Grho .* LargeScale .* abs(fading).^2;
                        for p = 1
                            sigPower_temp(ran,p) = power_temp(p,p);
                            for q = [1:p-1 p+1:num_UEs]
                                intfPower_temp(ran,p) = intfPower_temp(ran,p) + power_temp(p,q);
                            end
                        end
                    end
                elseif is_3Dbeamforming
                   for ran = 1:num_fading
                       for p = 1
                            Hi = ch_3D_temp(:,:,p,ran);
                            sigPower_temp(ran,p) = P_t * norm(Hi(:,p))^2; % MRT beamforming
                            for q = [1:p-1 p+1:num_UEs]
                                h_que = ch_3D_temp(:,q,p+1,ran);
                                w_que = h_que/norm(h_que); % MRT beamforming
                                h_inf = Hi(:,q);
                                intfPower_temp(ran,p) = intfPower_temp(ran,p) + P_t * abs(w_que'*h_inf)^2;
                            end
                       end
                   end
                end
                
                SINR = mean( sigPower_temp ./ (intfPower_temp + NoisePower) );
                barS_bench_temp(nt,np) = mean( sigPower_temp );
                barI_bench_temp(nt,np) = mean( intfPower_temp );
                R_temp_bench(nt,np) = log2(1+SINR)';
            end
        end
        barS(:,:,gi) = barS_irs_temp;
        barI(:,:,gi) = barI_irs_temp;
        barS(:,:,gi+2) = barS_bench_temp;
        barI(:,:,gi+2) = barI_bench_temp;
        
        mean_sigPow_lb(:,:,gi) = mean_sigPow_lb_temp;
        mean_itfPow_ub(:,:,gi) = mean_itfPow_ub_temp;
        
        R(:,:,gi) = R_temp;
        R(:,:,gi+2) = R_temp_bench;
    end
    
    barbarS(ue,:,1) = mean(barS(:,1,:),1);
    barbarS(ue,:,2) = mean(barS(:,2,:),1);
    
    barbarI(ue,:,1) = mean(barI(:,1,:),1);
    barbarI(ue,:,2) = mean(barI(:,2,:),1);
    
    avg_mean_SigPow_lb(ue,:,1) = mean(mean_sigPow_lb(:,1,:),1);
    avg_mean_SigPow_lb(ue,:,2) = mean(mean_sigPow_lb(:,2,:),1);
    
    avg_mean_ItrPow_ub(ue,:,1) = mean(mean_itfPow_ub(:,1,:),1);
    avg_mean_ItrPow_ub(ue,:,2) = mean(mean_itfPow_ub(:,2,:),1);
    
    avgThroughput(ue,:,1) = mean(R(:,1,:),1);
    avgThroughput(ue,:,2) = mean(R(:,2,:),1);
    
    toc
end

%% 5) save data
date = char(datetime('today', 'InputFormat','yyyy-MM-dd'));
save(['../data/multiCells/', date, '_barbarSinMultipleCells'],'barbarS');
save(['../data/multiCells/', date, '_barbarIinMultipleCells'],'barbarI');
save(['../data/multiCells/', date, '_AvgMeanSigPowlbinMultipleCells'],'avg_mean_SigPow_lb');
save(['../data/multiCells/', date, '_AvgMeanItrPowubinMultipleCells'],'avg_mean_ItrPow_ub');
save(['../data/multiCells/', date, '_avgThroughputinMultipleCells'],'avgThroughput');
