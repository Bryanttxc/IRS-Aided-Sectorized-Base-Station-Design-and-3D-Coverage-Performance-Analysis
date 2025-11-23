% created by Chen Xintong on 2022-5-3
% modified on 2023-12-17
% modified on 2025-11-22, mainly reFormat code
% -----------------------------------------------------------------------
% note: Please revise "theta_g" and "theta_m" in para_init.m first before
% running the code
% -----------------------------------------------------------------------

clear;
clc;
close all;

addpath("Geometry");
addpath("ERP");
addpath("Pathloss");

para_init;

%% 1) Generate locations of TX,IRS,UE and angles
num_IRS_idx = 3;
[UE_2D_raw, UE_2D, sel_idx, dis_irs_tx_3D, b, theta_t, theta_tx] = CalLocationAndAngle(num_IRS_idx);

%% 2) K-factor & Large-scale pathloss
num_UEs = size(UE_2D,1); % number of UEs
hue_idx = 1;
UE_z = Hue(hue_idx);
H_u = UE_z * ones(num_UEs,1);

% 2) -1] distances & angles
% Sector represents the IRS-aided BS including TX and IRS
% 2D, 3D distance from sector to each UE
sector_2D = [0 0];
vec_sector2UE_2D = UE_2D-sector_2D; %% can be negative
vec_sector2UE_z = H_u-H_i; %% can be negative
dis_sector2UE_2D = sqrt(vec_sector2UE_2D.^2);
dis_sector2UE_3D = sqrt(sum(dis_sector2UE_2D.^2,2)+vec_sector2UE_z.^2);

% IRS-UE angle, theta_r
gi = 1;
mat_sector2UE = [vec_sector2UE_2D vec_sector2UE_z];
theta_r = acosd( mat_sector2UE*b' ./ dis_sector2UE_3D );
F_r_LOS = GetRadiationValue(theta_r, G_i(gi));

% 2) -2] LoS Probability & LargeScale Pathloss
Prob_LOS = zeros(num_UEs,num_sectors); % LOS Probability
PL_overall = zeros(num_UEs,num_sectors); % Pathloss
is_LOS = zeros(num_UEs,num_sectors);
for p = 1:num_UEs
    for q = 1:num_sectors
        Prob_LOS(p,q) = LOSprobability(UE_z,dis_sector2UE_2D(p,q));
        if rand(1) < Prob_LOS(p,q) % LOS
            PL_overall(p,q) = LOSpathloss(UE_z,dis_sector2UE_2D(p,q),dis_sector2UE_3D(p,q),H_i,f_c,c);
            is_LOS(p,q) = 1;
        else % NLOS
            PL_overall(p,q) = NLOSpathloss(UE_z,dis_sector2UE_2D(p,q),dis_sector2UE_3D(p,q),H_i,f_c,c);
            is_LOS(p,q) = 0;
        end
    end
end
LargeScale_dB = -PL_overall;
LargeScale = 10.^(LargeScale_dB ./ 10);

% 2) -3] Rician K-factor
if UE_z <= H_i  % 3GPP TR25.996 formula for ground UE
    Kfactor = 10.^((13-0.03*dis_sector2UE_3D)./10);
else            %《3D trajectory ...》 for aerial UE
    theta_K = asin(abs(H_u-H_i)./dis_sector2UE_3D);
    Kfactor = funcRicianK(theta_K);
end
idx_NLoS = find(is_LOS == 0); % NLOS path
Kfactor(idx_NLoS) = 0;

% Generate multipath channel power
UE_idx = 1;
num_trial_angle = 10;
num_trial_phase = 300;
path = [3 30];
delta_cos = 1e-6;
delta_von = 1e-6;

Channel_MC     = zeros(num_trial_angle*num_trial_phase,length(path),length(theta_g));
Channelprim_MC = zeros(num_trial_angle*num_trial_phase,length(path),length(theta_g));
for ng = 1:length(theta_g)

    %% 3) Power gain of scattered waves
    [EA, ~] = CalPowerGainofScatteredWave(ng, num_IRS_idx);
    
    %% 4) Monte Carlo simulation
    for np = 1:length(path)
        
        num_path = path(np); 
        
        for tri_angle_idx = 1:num_trial_angle
            theta_MC = zeros(num_path,1);
            phi_MC = zeros(num_path,1);
            cosine_pdf_MC = zeros(num_path,1);
            vonMises_pdf_MC = zeros(num_path,1);
            cos_pdf_path_cnt = 0; 
            vonMise_pdf_path_cnt = 0;
            
            % Acceptance-Rejection Method
            % Step 1: generate a random number X uniformly distributed in the
            %         target domain
            % Step 2: calculate its corresponding target PDF value f(X)
            % Step 3: generate a random number U uniformly distributed in [0,m]
            %         where m = max(f(x))
            % Step 4: if U <= f(X), X will be chosen as the candidate number
            % Repeat Step 1-4 until satisfying the exit condition
            
            % vertical distribution 
            while cos_pdf_path_cnt<num_path
                temp_theta = (upper_cos(ng)-lower_cos(ng))*rand(1) + lower_cos(ng);
                temp_cosinePDF = funcCosinePdf(temp_theta, theta_g(ng), theta_m(ng));  
                randTheta = max_cosine_pdf(ng) * rand(1);
                if randTheta <= temp_cosinePDF
                    cos_pdf_path_cnt = cos_pdf_path_cnt + 1;
                    theta_MC(cos_pdf_path_cnt) = temp_theta;
                    cosine_pdf_MC(cos_pdf_path_cnt) = temp_cosinePDF;
                end
            end
            
            % horizontal distribution
            while vonMise_pdf_path_cnt<num_path
                temp_phi = (upper_von-lower_von)*rand(1) + lower_von;
                temp_vonMises_pdf = funcVonMisesPdf(temp_phi, phi_mu, k);
                rand_phi = max_vonMises_pdf * rand(1);
                if rand_phi <= temp_vonMises_pdf
                    vonMise_pdf_path_cnt = vonMise_pdf_path_cnt + 1;
                    phi_MC(vonMise_pdf_path_cnt) = temp_phi;
                    vonMises_pdf_MC(vonMise_pdf_path_cnt) = temp_vonMises_pdf;
                end
            end
    
            % % check the correction of random generation
            % num = 50;
            % [x,c] = hist(phi_MC,num);
            % dc = 2*pi/num;
            % x = x/numPath/dc;
            % bar(c,x,1)
            % tt = linspace(-pi,pi,1000);
            % ff = exp(k*cos(tt-phi_mu))/(2*pi*besseli(0,k));
            % hold on
            % plot(tt,ff)
            
            power_distrib_multipath = 1/(Kfactor(UE_idx)+1) * (cosine_pdf_MC * delta_cos) .* (vonMises_pdf_MC * delta_von);
            normal_pdm = power_distrib_multipath ./ sum(power_distrib_multipath); % normalized
            
            % radiation pattern
            theta_r_MC = theta_MC./pi*180;
            F_r_NLOS_theta_MC = GetRadiationValue(theta_r_MC, G_i(gi));
            
            h_LoS = exp(1j*2*pi.*dis_sector2UE_3D(UE_idx)./wavelength); % |h_LoS| = 1
            h_LoSprim = h_LoS .* sqrt(G_i(gi) .* F_r_LOS(UE_idx));
            
            for ntP = 1:num_trial_phase
                h_NLoS = sqrt(normal_pdm).*exp(1j*2*pi*rand(num_path,1));
                h_NLoSprim = h_NLoS .* sqrt(G_i(gi) .* F_r_NLOS_theta_MC); % multiply own coefficients
    
                Channel_MC((tri_angle_idx-1)*num_trial_phase+ntP,np,ng) =  sqrt(LargeScale(UE_idx)) * ( sqrt(Kfactor(UE_idx)./(Kfactor(UE_idx)+1)).*h_LoS + sqrt(1./(Kfactor(UE_idx)+1)).*sum(h_NLoS) );
                Channelprim_MC((tri_angle_idx-1)*num_trial_phase+ntP,np,ng) =  sqrt(LargeScale(UE_idx)) * ( sqrt(Kfactor(UE_idx)./(Kfactor(UE_idx)+1)).*h_LoSprim + sqrt(1./(Kfactor(UE_idx)+1)).*sum(h_NLoSprim) );
            end
        end
    end
    
    % theory curve
    idx_LoS = find(is_LOS == 1); % LOS path
    G_k = G_i(gi) .* F_r_LOS ./ EA(gi); % Rician K-factor gain
    
    KfactorPrim = zeros(num_UEs,num_sectors); % Kfactor considering radiation pattern
    G_rho = EA(gi).*ones(num_UEs,num_sectors); % SNR gain
    v2 = zeros(num_UEs,num_sectors);
    sigma2 = 1/2.*ones(num_UEs,num_sectors);
    
    KfactorPrim(idx_LoS) = G_k(idx_LoS) .* Kfactor(idx_LoS);
    G_rho(idx_LoS) = Kfactor(idx_LoS) ./ (Kfactor(idx_LoS)+1) .*G_i(gi).*F_r_LOS(idx_LoS) + 1./(Kfactor(idx_LoS)+1) .*EA(gi);
    v2(idx_LoS) = KfactorPrim(idx_LoS)./(KfactorPrim(idx_LoS)+1); % noncentral factor
    sigma2(idx_LoS) = 1./(2.*(KfactorPrim(idx_LoS)+1));
    
    % method 1
    x_temp{ng} = min(min(abs(Channelprim_MC(:,:,ng)))):1e-7:max(max(abs(Channelprim_MC(:,:,ng))));
    x = x_temp{ng} ./ sqrt(LargeScale(UE_idx)*G_rho(UE_idx));
    Channelprim_Theory{ng} = cdf('Rician', x, sqrt(v2(UE_idx)), sqrt(sigma2(UE_idx)));
    % method 2
    % a = sqrt(v2(UEidx)) / sqrt(sigma2(UEidx));
    % b = x / sqrt(sigma2(UEidx));
    % Channelprim_Theory = 1 - marcumq(a,b,1);

end

%% 5) save data
date = char(datetime('today', 'InputFormat','yyyy-MM-dd'));
save(['../data/ricianProve/', date, '_channelPrim'],"Channelprim_MC", "x_temp", "Channelprim_Theory");

%% 6) plot
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

figure;
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

% Theory curve without ERP
% temp = min(min(abs(Channel_MC))):1e-8:max(max(abs(Channel_MC)));
% x = temp ./ sqrt(LargeScale(UEidx));
% v22 = Kfactor(UEidx) ./ (Kfactor(UEidx)+1);
% sigma22 = 1 ./ (2*(Kfactor(UEidx)+1));
% Channel_theory = cdf('Rician', x, sqrt(v22), sqrt(sigma22));
% h_NLoS_theory = cdf('Rayleigh',temp,sqrt(sum(NormalPDM)/2));

% plot without ERP
% h7 = plot(temp{1},Channel_Theory{1});
% set(h7,'linestyle','-','linewidth',linewidth_theory,'color',[1.00,0.32,0.41]) % [0.01,0.62,0.65]
% hold on
% h8 = cdfplot(10*log10(abs(Channel_MC(:,1,1)).^2));
% h8.XData = h8.XData(index2);
% h8.YData = h8.YData(index2);
% set(h8,'linestyle',':','linewidth',linewidth_MC,'Marker','^','MarkerSize',MarkerSize_MC,'MarkerFaceColor',color(6,:),'color',color(6,:)) % [0.00,0.45,0.74]
% hold on
% h9 = cdfplot(10*log10(abs(Channel_MC(:,2,1)).^2));
% h9.XData = h9.XData(index2);
% h9.YData = h9.YData(index2);
% set(h9,'linestyle','-.','linewidth',linewidth_MC,'Marker','s','MarkerSize',MarkerSize_MC,'MarkerFaceColor','k','color','k')
% hold on
% 
% h10 = plot(temp{2},Channel_Theory{2});
% set(h10,'linestyle','-','linewidth',linewidth_theory,'color',color(3,:));
% hold on
% h11 = cdfplot(10*log10(abs(Channel_MC(:,1,2)).^2));
% h11.XData = h11.XData(index2);
% h11.YData = h11.YData(index2);
% set(h11,'linestyle',':','linewidth',linewidth_MC,'Marker','^','MarkerSize',MarkerSize_MC,'MarkerFaceColor','g','color','g'); % [0.00,0.45,0.74]
% hold on
% h12 = cdfplot(10*log10(abs(Channel_MC(:,2,2)).^2));
% h12.XData = h12.XData(index2);
% h12.YData = h12.YData(index2);
% set(h12,'linestyle','-.','linewidth',linewidth_MC,'Marker','s','MarkerSize',MarkerSize_MC,'MarkerFaceColor',[0.89,0.35,1.00],'color',[0.89,0.35,1.00]);%[0.01,0.62,0.65]
