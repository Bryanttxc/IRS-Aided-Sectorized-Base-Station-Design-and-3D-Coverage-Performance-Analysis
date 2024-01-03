%created by Chen Xintong on 2023-12-17

clc;
clear;
close all

para_init;

%% Generate sector locations in multiple cells
[SectorLocation, LabelLocation, NumSectors] = GenerateSectorLoc();

%% Generate location of TX,IRS and UE & angle
nnn = 3;
[UE_2D_raw, UE_2D, sel_idx, dis_irs_tx_3D, b, theta_t, theta_tx] = calLocationAndAngle(nnn);

% boresight of each sector
% b1_i : sector 1 4 7 10 13 16 19
% b2_i : sector 2 5 8 11 14 17 20
% b3_i : sector 3 6 9 12 15 18 21

b1_i = [cosd(0)*cosd(alpha_y)   sind(0)*cosd(alpha_y)   sind(alpha_y)];
b2_i = [cosd(120)*cosd(alpha_y) sind(120)*cosd(alpha_y) sind(alpha_y)];
b3_i = [cosd(240)*cosd(alpha_y) sind(240)*cosd(alpha_y) sind(alpha_y)];

%% Power gain of scattered waves
[EA, EC] = calPowerGainofScatteredWave(1,nnn);

%% Ergodic Throughput %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumRow = 2; % Row1: the channel between the UE in cell 1 and all sectors, Row2: the channel between each sector and its own served UE 
NumUEs = NumSectors; % assume each Sector have one UE
R_temp = zeros(numTrials,numPlot); % average Throughput
R_temp_bench = zeros(numTrials,numPlot); % average Throughput
R = zeros(numTrials,numPlot,numIRSscheme+numBenchmark); % Four scheme
avgThroughput = zeros(size(UE_2D,1),numIRSscheme+numBenchmark,numPlot);

barS_irs_temp = zeros(numTrials,numPlot);
barI_irs_temp = zeros(numTrials,numPlot);
barS_bench_temp = zeros(numTrials,numPlot);
barI_bench_temp = zeros(numTrials,numPlot);

barS = zeros(numTrials,numPlot,numIRSscheme+numBenchmark); % Four scheme
barI = zeros(numTrials,numPlot,numIRSscheme+numBenchmark);

barbarS = zeros(size(UE_2D,1),numIRSscheme+numBenchmark,numPlot);
barbarI = zeros(size(UE_2D,1),numIRSscheme+numBenchmark,numPlot);

MeanSigPowlb_temp = zeros(numTrials,numPlot); % lower bound
MeanItrPowub_temp = zeros(numTrials,numPlot); % upper bound
MeanSigPowlb = zeros(numTrials,numPlot,numIRSscheme+numBenchmark); % Four scheme
MeanItrPowub = zeros(numTrials,numPlot,numIRSscheme+numBenchmark);
AvgMeanSigPowlb = zeros(size(UE_2D,1),numIRSscheme+numBenchmark,numPlot);
AvgMeanItrPowub = zeros(size(UE_2D,1),numIRSscheme+numBenchmark,numPlot);

% need very long time
for num = 1:size(UE_2D,1)
    num
    
    tic
    %%%%%%%%%%%%%%% Generate UE random locations of other Sectors %%%%%%%%%%%%%%%
    UEmaxRadius = r * sqrt(3)/2; % inscribed circle
    theta = rand(NumUEs-1,numTrials)* 2*pi; % uniform distribution within the circle
    r_ue = rand(NumUEs-1,numTrials).* UEmaxRadius;
    x_ue = r_ue.*cos(theta) + LabelLocation(2:end,1);
    y_ue = r_ue.*sin(theta) + LabelLocation(2:end,2);

    OtherSector_UE_Point(:,:,1) = x_ue;
    OtherSector_UE_Point(:,:,2) = y_ue;
    
    Sel_UE_Point(:,:,1) = [repmat(UE_2D(num,1),1,numTrials);OtherSector_UE_Point(:,:,1)];
    Sel_UE_Point(:,:,2) = [repmat(UE_2D(num,2),1,numTrials);OtherSector_UE_Point(:,:,2)];
    
    for gi = 1:length(G_i)
        G_i_temp = G_i(gi);
        EA_temp = EA(gi);
        EC_temp = EC(gi);
        
        %%%%%% radiation pattern %%%%%%
        F_tx = getRadiationValue(theta_tx, G_tx);
        F_t  = getRadiationValue(theta_t, G_i_temp);
        [MinFcombine,Minidx] = min(F_tx.*F_t);
        [MaxFcombine,Maxidx] = max(F_tx.*F_t);

        parfor nt = 1:numTrials

            %%% 2D,3D distance and angleH, angleV from each UE to each sector %%%
            vec_SectorToUE_x = bsxfun(@minus, Sel_UE_Point(:,nt,1), SectorLocation(:,1)'); %% can be negative
            vec_SectorToUE_y = bsxfun(@minus, Sel_UE_Point(:,nt,2), SectorLocation(:,2)'); %% can be negative

            %%% delete unnecessary data %%%
            Sel_vec_SectorToUE_x = [vec_SectorToUE_x(1,:);diag(vec_SectorToUE_x)'];
            Sel_vec_SectorToUE_y = [vec_SectorToUE_y(1,:);diag(vec_SectorToUE_y)'];

            DistToSector2D = sqrt(Sel_vec_SectorToUE_x.^2+Sel_vec_SectorToUE_y.^2);

            %%% angleH of Fixed BS Pattern & 3D beamforming %%%
            angleHatSectorGlobal = atan2d(Sel_vec_SectorToUE_y,Sel_vec_SectorToUE_x); % degree
            angleHatSectorLocal = angleHatSectorGlobal; %% relative to the boresight direction of each sector, 0,120,-120
            angleHatSectorLocal(:,2:3:end) = angleHatSectorGlobal(:,2:3:end)-120;
            angleHatSectorLocal(:,3:3:end) = angleHatSectorGlobal(:,3:3:end)-(-120);

            %%% convert to -180~180 %%%
            idxNegt = find(angleHatSectorLocal<-180);
            angleHatSectorLocal(idxNegt) = angleHatSectorLocal(idxNegt)+360;
            idxGreater180 = find(angleHatSectorLocal>180);
            angleHatSectorLocal(idxGreater180) = angleHatSectorLocal(idxGreater180)-360;
            index = find(angleHatSectorLocal < -90 | angleHatSectorLocal > 90);

            for np = 1:numPlot

                UE_z = Hue(np)*ones(NumUEs,1);
                vec_SectorToUE_z = bsxfun(@minus, UE_z, H_i*ones(NumSectors,1)'); %% can be negative
                Sel_vec_SectorToUE_z = [vec_SectorToUE_z(1,:);diag(vec_SectorToUE_z)'];
                DistToSector3D = sqrt(DistToSector2D.^2+Sel_vec_SectorToUE_z.^2);

                %%% The angleV of FixedPattern & 3D beamforming BS %%%
                angleVatSectorLocal=acosd(Sel_vec_SectorToUE_z./DistToSector3D);

                %%%%%%%%%%%%%%%%%%%% The angleV of cos & cos3 IRS %%%%%%%%%%%%%%%%%%%%
                Sector1ToUE = [reshape(Sel_vec_SectorToUE_x(:,1:3:end),[],1) reshape(Sel_vec_SectorToUE_y(:,1:3:end),[],1) reshape(Sel_vec_SectorToUE_z(:,1:3:end),[],1)];
                Sector2ToUE = [reshape(Sel_vec_SectorToUE_x(:,2:3:end),[],1) reshape(Sel_vec_SectorToUE_y(:,2:3:end),[],1) reshape(Sel_vec_SectorToUE_z(:,2:3:end),[],1)]; 
                Sector3ToUE = [reshape(Sel_vec_SectorToUE_x(:,3:3:end),[],1) reshape(Sel_vec_SectorToUE_y(:,3:3:end),[],1) reshape(Sel_vec_SectorToUE_z(:,3:3:end),[],1)];

                theta_r = zeros(NumRow,NumSectors);
                theta_r(:,1:3:end) = reshape( acosd( Sector1ToUE*b1_i' ./ reshape(DistToSector3D(:,1:3:end),[],1) )' ,NumRow,[]);
                theta_r(:,2:3:end) = reshape( acosd( Sector2ToUE*b2_i' ./ reshape(DistToSector3D(:,2:3:end),[],1) )' ,NumRow,[]);
                theta_r(:,3:3:end) = reshape( acosd( Sector3ToUE*b3_i' ./ reshape(DistToSector3D(:,3:3:end),[],1) )' ,NumRow,[]);

                %%%%%%%%%%%%%%%  LOS Probability and Large-Scale Pathloss %%%%%%%%%%%%%%%
                PrLOS = zeros(NumRow,NumSectors); % LOS Probability
                PL_overall = zeros(NumRow,NumSectors);
                isLOS = zeros(NumRow,NumSectors);
                for p = 1:NumRow
                    for q = 1:NumSectors
                        PrLOS(p,q) = LOSprobability(UE_z(1),DistToSector2D(p,q));
                    end
                end
                PrLOS(index) = 0; % modify
                
                for p = 1:NumRow
                    for q = 1:NumSectors
                        if rand(1) < PrLOS(p,q) %% LOS
                            PL_overall(p,q) = LOSpathloss(UE_z(1),DistToSector2D(p,q),DistToSector3D(p,q),H_i,f_c,c);
                            isLOS(p,q) = 1;
                        else %% NLOS
                            PL_overall(p,q) = NLOSpathloss(UE_z(1),DistToSector2D(p,q),DistToSector3D(p,q),H_i,f_c,c);
                            isLOS(p,q) = 0;
                        end
                    end
                end
                LargeScaledB = -PL_overall;
                LargeScale = 10 .^(LargeScaledB ./ 10);

                %%%%%%%%%%%%%%% Rician Kfactor %%%%%%%%%%%%%%%
                if UE_z(1) <= H_i
                    Kfactor_LOS = 10.^((13-0.03*DistToSector3D)./10);
                else
                    theta_K = asin(abs(UE_z(1)'-H_i)./DistToSector3D);
                    Kfactor_LOS = funcRicianK(theta_K);
                end
                Kfactor_LOS(index) = 0;
                Kfactor_NLOS = zeros(NumRow,NumSectors); 

                %%%%%%%%%%%%%%%  IRS cos & cos3 pattern %%%%%%%%%%%%%%%
                F_r = getRadiationValue(theta_r, G_i_temp);

                % Rician K-factor gain
                Gk = G_i_temp .* F_r ./ EA_temp; 

                %%%%%% determine parameters of fading %%%%%%
                KfactorPrim = zeros(NumRow,NumSectors); % Kfactor considering antenna radiation pattern
                Grho = zeros(NumRow,NumSectors); % SNR gain
                a2 = zeros(NumRow,NumSectors);
                sigma2 = zeros(NumRow,NumSectors);
                for p = 1
                    for q = 1:NumSectors
                        if isLOS(p,q) %% LOS
                            Grho(p,q) = Kfactor_LOS(p,q)./(Kfactor_LOS(p,q)+1) .*G_i(gi).*F_r(p,q) + 1./(Kfactor_LOS(p,q)+1) .*EA_temp;
                            KfactorPrim(p,q) = Gk(p,q) .* Kfactor_LOS(p,q);
                            a2(p,q) = KfactorPrim(p,q)./( 2.*(KfactorPrim(p,q)+1) ); % noncentral factor
                            sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                        else %%NLOS
                            Grho(p,q) = EA_temp;
                            KfactorPrim(p,q) = Kfactor_NLOS(p,q);
                            a2(p,q) = KfactorPrim(p,q)./( 2.*(KfactorPrim(p,q)+1) ); %noncentral factor
                            sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                        end
                    end
                end

                %%%%%%% Generate random channel & SINR %%%%%%%
                SignalPower_temp = zeros(numfading,1);
                IntrfrPower_temp = zeros(numfading,1);
                for ran = 1:numfading
                    for p = 1
                        fading = ( sqrt(a2(p,p))+sqrt(sigma2(p,p)).*randn(M(nnn)*N(nnn),1) )+...
                            1j*( sqrt(a2(p,p))+sqrt(sigma2(p,p)).*randn(M(nnn)*N(nnn),1) );% small scale fading
                        SignalPower_temp(ran,p) =  P_t * G_tx*G_i(gi)*beta0*A^2.*Grho(p,p) .* LargeScale(p,p).*...
                                      abs( sum( abs(fading).* sqrt(F_tx.*F_t) ./ dis_irs_tx_3D ) ).^2;
                        for q = [1:p-1 p+1:NumUEs] % for other UEs
                            randPhase = exp(1j*2*pi*rand(M(nnn)*N(nnn),1));
                            fading = ( sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M(nnn)*N(nnn),1) )+...
                                1j*( sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M(nnn)*N(nnn),1) );% small scale fading
                            IntrfrPower_temp(ran,p) = IntrfrPower_temp(ran,p) + P_t * G_tx*G_i(gi)*beta0.*Grho(p,q).* LargeScale(p,q).*...
                                  abs( sum( abs(fading).* sqrt(F_tx.*F_t).*A.*randPhase.*exp(-1j*2*pi/wavelength.*dis_irs_tx_3D).*exp(-1j*2*pi/wavelength.*DistToSector3D(p,q))./dis_irs_tx_3D ) ).^2;
                        end
                    end
                end
                SINR = mean( SignalPower_temp ./ (IntrfrPower_temp + NoisePower) );
                barS_irs_temp(nt,np) = mean(SignalPower_temp);
                barI_irs_temp(nt,np) = mean(IntrfrPower_temp);
                MeanSigPowlb_temp(nt,np) = ( pi / (4*(KfactorPrim(1,1)+1))*L(-KfactorPrim(1,1))^2*(M(nnn)*N(nnn))^2 + (1 - pi/(4*(KfactorPrim(1,1)+1)))*L(-KfactorPrim(1,1))^2*M(nnn)*N(nnn) ).*...
                                        P_t* (A^2*G_tx*G_i(gi)*Grho(1,1)*MinFcombine*beta0/dis_irs_tx_3D(Minidx)*LargeScale(1,1));
                MeanItrPowub_temp(nt,np) = sum( M(nnn)*N(nnn)* P_t* (A^2*G_tx*G_i(gi)*MaxFcombine*beta0/dis_irs_tx_3D(Maxidx).*Grho(1,2:end).*LargeScale(1,2:end)) );
                R_temp(nt,np) = log2(1+SINR)';


               %% Fixed BS Pattern & 3D beamforming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    isFixedPattern = 0;
                    is3Dbeamforming = 0;

                    %%%%%% Sector antennas' Gain %%%%%%
                    if gi+2 == 3 %% FixedPattern
                        isFixedPattern = 1;
                        GsectordB = zeros(NumRow,NumSectors);
                        for p = 1:NumRow
                            for q = 1:NumSectors
                                GsectordB(p,q) = ArrayPowerPattern(angleVatSectorLocal(p,q),angleHatSectorLocal(p,q),angleTiltV,angleTiltH,M_b(nnn),N_b(nnn),dV,dH,wavelength);
                            end
                        end
                    elseif gi+2 == 4 %% 3D beamforming
                        is3Dbeamforming = 1;
                        GsectordB = zeros(NumRow,NumSectors);
                        for p = 1:NumRow
                            for q = 1:NumSectors
                                GelementdB = ElementPowerPatternOverall(angleVatSectorLocal(p,q),angleHatSectorLocal(p,q)); % the element power gain
                                Gelement = 10^(GelementdB/10);
                                Felement = sqrt(Gelement); % the element field pattern
                                rMatrix = zeros(M_b(nnn)*N_b(nnn),3);
                                idx = -N_b(nnn)/2+1:N_b(nnn)/2;
                                cnt = 1;
                                for v = 1:N_b(nnn)
                                     rMatrix(1+(cnt-1)*M_b(nnn):cnt*M_b(nnn),3) = [-M_b(nnn)/2:(M_b(nnn)/2-1)]'*dV;
                                     rMatrix(1+(cnt-1)*M_b(nnn):cnt*M_b(nnn),2) = (idx(v)-1)*dH;
                                     cnt = cnt + 1;
                                end
                                kVector = (2*pi/wavelength)*[sind(angleVatSectorLocal(p,q))*cosd(angleHatSectorLocal(p,q)),sind(angleVatSectorLocal(p,q))*sind(angleHatSectorLocal(p,q)),cosd(angleVatSectorLocal(p,q))];%wave vector
                                SteeringVector = exp(-1i.*(rMatrix*kVector'));%Num_Ant x 1
                                GsectordB(p,q) = 10*log10( abs(Felement*SteeringVector(centeridx(nnn)))^2 );
                            end
                        end
                    end

                    Gk = 10.^(GsectordB./10) ./ EC_temp; % Ricial K-factor gain
                    Gk(index) = 0;

                    %%%%%% determine parameters of each UE %%%%%% 
                    Grho        = zeros(NumRow,NumSectors);
                    KfactorPrim = zeros(NumRow,NumSectors);
                    a2          = zeros(NumRow,NumSectors);
                    sigma2      = zeros(NumRow,NumSectors);
                    fading      = zeros(NumRow,NumSectors);
                    ch_3D_temp  = zeros(M_b(nnn)*N_b(nnn),NumSectors,NumRow,numfading);                
                    for p = 1:NumRow
                        for q = 1:NumSectors
                            if isLOS(p,q)
                                if isFixedPattern
                                    Grho(p,q) = ( Kfactor_LOS(p,q)./(Kfactor_LOS(p,q)+1) .*(10.^(GsectordB(p,q)./10)) + 1./(Kfactor_LOS(p,q)+1) .*EC_temp ) ./ (M_b(nnn)*N_b(nnn)); % SNR gain
                                elseif is3Dbeamforming
                                    Grho(p,q) = Kfactor_LOS(p,q)./(Kfactor_LOS(p,q)+1) .*(10.^(GsectordB(p,q)./10)) + 1./(Kfactor_LOS(p,q)+1) .*EC_temp; % SNR gain
                                end
                                KfactorPrim(p,q) = Gk(p,q) .* Kfactor_LOS(p,q);
                                a2(p,q) = KfactorPrim(p,q)./(2.*(KfactorPrim(p,q)+1)); % noncentral factor
                                sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                            else
                                if isFixedPattern
                                    Grho(p,q) = EC_temp / (M_b(nnn)*N_b(nnn)); % SNR gain
                                elseif is3Dbeamforming
                                    Grho(p,q) = EC_temp;
                                end
                                KfactorPrim(p,q) = Gk(p,q) .* Kfactor_NLOS(p,q);                            
                                a2(p,q) = KfactorPrim(p,q)./(2.*(KfactorPrim(p,q)+1)); % noncentral factor
                                sigma2(p,q) = 1./(2.*(KfactorPrim(p,q)+1));
                            end

                            if is3Dbeamforming
                                for ran = 1:numfading
                                     ch_3D_temp(:,q,p,ran) = sqrt(Grho(p,q)) .* sqrt(LargeScale(p,q)) .* ...
                                         ( (sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M_b(nnn)*N_b(nnn),1)) + 1j* (sqrt(a2(p,q))+sqrt(sigma2(p,q)).*randn(M_b(nnn)*N_b(nnn),1)) );
                                end
                            end

                        end
                    end

                    %%%%%% calculate each UE's SINR %%%%%%
                    PowerTemp = zeros(NumRow,NumSectors);
                    SignalPower_temp = zeros(numfading,1);
                    IntrfrPower_temp = zeros(numfading,1);
                    if isFixedPattern
                        for ran = 1:numfading
                            fading = ( sqrt(a2)+sqrt(sigma2).*randn(NumRow,NumSectors) ) +...
                                        1j* ( sqrt(a2)+sqrt(sigma2).*randn(NumRow,NumSectors) );
                            PowerTemp = P_t .* Grho .* LargeScale .* abs(fading).^2;
                            for p = 1
                                SignalPower_temp(ran,p) = PowerTemp(p,p);
                                for q = [1:p-1 p+1:NumUEs]
                                    IntrfrPower_temp(ran,p) = IntrfrPower_temp(ran,p) + PowerTemp(p,q);
                                end
                            end
                        end
                    elseif is3Dbeamforming
                       for ran = 1:numfading
                           for p = 1
                                Hi = ch_3D_temp(:,:,p,ran);
                                SignalPower_temp(ran,p) = P_t * norm(Hi(:,p))^2; %% MRT beamforming
                                for q = [1:p-1 p+1:NumUEs]
                                    h_que = ch_3D_temp(:,q,p+1,ran);
                                    w_que = h_que/norm(h_que); %% MRT beamforming
                                    h_inf = Hi(:,q);
                                    IntrfrPower_temp(ran,p) = IntrfrPower_temp(ran,p) + P_t * abs(w_que'*h_inf)^2;
                                end
                           end
                       end
                    end
                    
                    SINR = mean( SignalPower_temp ./ (IntrfrPower_temp + NoisePower) );
                    barS_bench_temp(nt,np) = mean( SignalPower_temp );
                    barI_bench_temp(nt,np) = mean( IntrfrPower_temp );
                    R_temp_bench(nt,np) = log2(1+SINR)';
               
            end
        end
        barS(:,:,gi) = barS_irs_temp;
        barI(:,:,gi) = barI_irs_temp;
        barS(:,:,gi+2) = barS_bench_temp;
        barI(:,:,gi+2) = barI_bench_temp;
        
        MeanSigPowlb(:,:,gi) = MeanSigPowlb_temp;
        MeanItrPowub(:,:,gi) = MeanItrPowub_temp;
        
        R(:,:,gi) = R_temp;
        R(:,:,gi+2) = R_temp_bench;
    end
    
    barbarS(num,:,1) = mean(barS(:,1,:),1);
    barbarS(num,:,2) = mean(barS(:,2,:),1);
    
    barbarI(num,:,1) = mean(barI(:,1,:),1);
    barbarI(num,:,2) = mean(barI(:,2,:),1);
    
    AvgMeanSigPowlb(num,:,1) = mean(MeanSigPowlb(:,1,:),1);
    AvgMeanSigPowlb(num,:,2) = mean(MeanSigPowlb(:,2,:),1);
    
    AvgMeanItrPowub(num,:,1) = mean(MeanItrPowub(:,1,:),1);
    AvgMeanItrPowub(num,:,2) = mean(MeanItrPowub(:,2,:),1);
    
    avgThroughput(num,:,1) = mean(R(:,1,:),1);
    avgThroughput(num,:,2) = mean(R(:,2,:),1);
    
    toc
end

%% save
save('../data/multiCells/20231217barbarSinMultipleCells','barbarS');
save('../data/multiCells/20231217barbarIinMultipleCells','barbarI');
save('../data/multiCells/20231217AvgMeanSigPowlbinMultipleCells','AvgMeanSigPowlb');
save('../data/multiCells/20231217AvgMeanItrPowubinMultipleCells','AvgMeanItrPowub');
save('../data/multiCells/20231217avgThroughputinMultipleCells','avgThroughput');
