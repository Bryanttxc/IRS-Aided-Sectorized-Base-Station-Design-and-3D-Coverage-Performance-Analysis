function [EA, EC] = calPowerGainofScatteredWave(gidx, numidx)
% CALPOWERGAINOFSCATTEREDWAVE calculate the power gain of scattered wave
% considering the radiation pattern effect
% Parameter explain:
% 1. gidx: index of mean scattered angles theta_g
% 2. numidx: index of the IRS element number
% 3. EA: power gain of scattered wave in IRS cases
% 4. EC: power gain of scattered wave in traditional antennas cases

para_init;

if ~isnan(numidx)
    M_b = M_b(numidx);
    N_b = N_b(numidx);
end

if ~isnan(gidx) % For RicianKfactorProve.m
    theta_m = theta_m(gidx);
    theta_g = theta_g(gidx);
    lower_cos = lower_cos(gidx);
    upper_cos = upper_cos(gidx);
end

% angular spectrum
% IRS-aided BS scheme
cosinePDF = pi/(4*theta_m) * cos(pi/(2*theta_m)*(pi/2-theta-theta_g)); % vertical distribution
vonMisesPDF = exp(k*cos(phi-phi_mu))/(2*pi*besseli(0,k)); % truncated horizontal distribution
normalVonMisesPDF = vonMisesPDF / eval(vpa(int(vonMisesPDF, phi, lower_von, upper_von))); % normalized

% Benchmark
cosinePDF_bs = 45/theta_m_bs.*cosd(90/theta_m_bs.*(90-theta_bs_ori-theta_g_bs));
normalCosinePDF_bs = cosinePDF_bs / sum(cosinePDF_bs);
normalCosinePDF_bs = repmat(normalCosinePDF_bs, length(phi_bs_ori), 1);

vonMisesPDF_bs = exp(k.*cosd(phi_bs_ori-phi_mu_bs)) ./ (2*pi*besseli(0,k));
normalVonMisesPDF_bs = vonMisesPDF_bs / sum(vonMisesPDF_bs); % normalized, sum = 1
normalVonMisesPDF_bs = reshape(repmat(normalVonMisesPDF_bs', length(theta_bs_ori), 1), [], 1 );

% impact analysis of antenna pattern
EA = zeros(numIRSscheme,1); % power gain of scattered waves for IRS-aided BS
EC = zeros(numBenchmark,1); % power gain of scattered waves for Benchmark
for n = 1: numIRSscheme+numBenchmark

    if n <= numIRSscheme % IRS-aided BS scheme
        F = cos(theta)^(G_i(n)/2-1); % ERP
        Ps_ip = eval(int(cosinePDF,theta,lower_cos,upper_cos)).* eval(vpa(int(normalVonMisesPDF,phi,lower_von,upper_von),4)); % isotropic pattern
        Ps_nip = eval(int(G_i(n)*F*cosinePDF,theta,lower_cos,upper_cos)).* eval(vpa(int(normalVonMisesPDF,phi,lower_von,upper_von),4)); % nonisotropic pattern
        EA(n) = Ps_nip ./ Ps_ip;

    else % benchmark schemes
        if n == numIRSscheme+1  % fixed Pattern
            ScatterPowerPatterndB = ArrayPowerPattern(theta_bs,phi_bs,angleTiltV,angleTiltH,M_b,N_b,dV,dH,wavelength);
        elseif n == numIRSscheme+2  % 3D beamforming
            ScatterPowerPatterndB = ElementPowerPatternOverall(theta_bs,phi_bs);
        end
        dS = 4*pi*((upper_cos_bs-lower_cos_bs)/360) / length(theta_bs); % dS = sin(theta)d(theta)d(phi)
        dThetadPhi = dS ./ sind(theta_bs-90);
        Ps_ip_bs = normalCosinePDF_bs .* normalVonMisesPDF_bs .* dThetadPhi;
        Ps_nip_bs = normalCosinePDF_bs .* normalVonMisesPDF_bs .* 10.^(ScatterPowerPatterndB./10) .* dThetadPhi;  
        EC(n-2) = sum(Ps_nip_bs(:)) / sum(Ps_ip_bs(:));

    end
    
end

