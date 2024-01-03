function [D_opt] = calDistanceFromTXtoIRS(N_i, hpbw, S_n)
% CALDISTANCEFROMTXTOIRS find the proper distance D to illuminate the
% whole IRS panel, TX antenna boresight is perpendicular to the IRS surface
% Parameter explain: 
% 1. N_i: number of IRS elements
% 2. hpbw: half-power beamwidth of TX beam
% 3. S_n: area of an IRS element

% D = [0.001:0.001:1]; % distance array
D = [0.001:0.001:5]; % distance array
theta_i = 0; % incident angle with respect to the IRS center
alpha = D .* sind(hpbw/2) / cosd(theta_i+hpbw/2); % semimajor axis a
gamma = sind(theta_i) / cosd(hpbw/2); % Eccentricity c/a
beta = alpha .* sqrt(1-gamma^2); % semiminor axis b
S = pi .* alpha .* beta; % IRS illuminated area
N_irs = floor(S ./ S_n);

D_opt = min( D(N_irs == N_i) );

end