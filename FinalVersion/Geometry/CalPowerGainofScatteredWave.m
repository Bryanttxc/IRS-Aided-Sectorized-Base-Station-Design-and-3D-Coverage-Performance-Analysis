function [EA, EC] = CalPowerGainofScatteredWave(g_idx, num_idx)
% CALPOWERGAINOFSCATTEREDWAVE calculate the power gain of scattered wave
% considering the radiation pattern effect
% Parameter explain:
% 1. gidx: index of mean scattered angles theta_g
% 2. numidx: index of the IRS element number
% 3. EA: power gain of scattered wave in IRS cases
% 4. EC: power gain of scattered wave in traditional antennas cases

para_init;

% cited by《A 3D geometry-based stochastic channel model for UAV-MIMO channels》
hdfunc_cosine_pdf = @(theta,theta_g,theta_m)pi/(4*theta_m) * cos(pi/(2*theta_m)*(pi/2-theta-theta_g));
hdfunc_vonMises_pdf = @(phi,phi_mu,k)exp(k*cos(phi-phi_mu))/(2*pi*besseli(0,k));

if ~isnan(num_idx)
    M_b = M_b(num_idx);
    N_b = N_b(num_idx);
end

if ~isnan(g_idx) % For RicianKfactorProve.m
    theta_m = theta_m(g_idx);
    theta_g = theta_g(g_idx);
    lower_cos = lower_cos(g_idx);
    upper_cos = upper_cos(g_idx);
end

%% 1) angular spectrum
% 1) -1] IRS-aided BS scheme, expression
cosine_pdf = hdfunc_cosine_pdf(theta,theta_g,theta_m); % vertical distribution
vonMises_pdf = hdfunc_vonMises_pdf(phi,phi_mu,k); % truncated horizontal distribution
normal_vonMises_pdf = vonMises_pdf / eval(vpa(int(vonMises_pdf, phi, lower_von, upper_von))); % normalized

% 1) -2] Benchmark, integrated
cosine_pdf_bs = hdfunc_cosine_pdf(theta_bs_ori_rad,theta_g_bs_rad,theta_m_bs_rad);
normal_cosine_pdf_bs = cosine_pdf_bs / sum(cosine_pdf_bs); % normalized
normal_cosine_pdf_bs = repmat(normal_cosine_pdf_bs, length(phi_bs_ori), 1);

normal_vonMises_pdf_bs = hdfunc_vonMises_pdf(phi_bs_ori_rad, phi_mu_bs_rad, k);
normal_vonMises_pdf_bs = normal_vonMises_pdf_bs / sum(normal_vonMises_pdf_bs); % normalized
normal_vonMises_pdf_bs = reshape(repmat(normal_vonMises_pdf_bs', length(theta_bs_ori), 1), [], 1);

%% 2) impact analysis of antenna pattern
EA = zeros(num_IRSscheme,1); % power gain of scattered waves for IRS-aided BS
EC = zeros(num_benchmark,1); % power gain of scattered waves for Benchmark
for scheme_idx = 1:num_IRSscheme+num_benchmark

    if scheme_idx <= num_IRSscheme % IRS-aided BS scheme
        F = cos(theta)^(G_i(scheme_idx)/2-1); % ERP
        Ps_ip = eval(int(cosine_pdf,theta,lower_cos,upper_cos)).* eval(vpa(int(normal_vonMises_pdf,phi,lower_von,upper_von),4)); % isotropic pattern
        Ps_nip = eval(int(G_i(scheme_idx)*F*cosine_pdf,theta,lower_cos,upper_cos)).* eval(vpa(int(normal_vonMises_pdf,phi,lower_von,upper_von),4)); % nonisotropic pattern
        EA(scheme_idx) = Ps_nip ./ Ps_ip;

    else % benchmark schemes
        if scheme_idx == num_IRSscheme+1  % fixed Pattern
            scatter_power_pattern_dB = ArrayPowerPattern(theta_bs,phi_bs,angleTiltV,angleTiltH,M_b,N_b,dV,dH,wavelength);
        elseif scheme_idx == num_IRSscheme+2  % 3D beamforming
            scatter_power_pattern_dB = ElementPowerPatternOverall(theta_bs,phi_bs);
        end
        dS = 4*pi*((upper_cos_bs-lower_cos_bs)/360) / length(theta_bs); % dS = sin(theta)d(theta)d(phi)
        dTheta_dPhi = dS ./ sind(theta_bs-90);
        Ps_ip_bs = normal_cosine_pdf_bs .* normal_vonMises_pdf_bs .* dTheta_dPhi;
        Ps_nip_bs = normal_cosine_pdf_bs .* normal_vonMises_pdf_bs .* 10.^(scatter_power_pattern_dB./10) .* dTheta_dPhi;  
        EC(scheme_idx-2) = sum(Ps_nip_bs(:)) / sum(Ps_ip_bs(:));
    end
    
end
