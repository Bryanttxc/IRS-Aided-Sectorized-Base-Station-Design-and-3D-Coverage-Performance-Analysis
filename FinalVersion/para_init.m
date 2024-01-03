%created by Chen Xintong on 2022-5-3

% General parameters
f_c           = 2e9;                                                % 2GHz Frequency
c             = 3e8;                                                % light speed
wavelength    = c/f_c;                                              % (m)
N0            = -174;                                               % NoiseDensity (dBm)
NF            = 6;                                                  % NoiseFigure (dB)
B             = 180e3;                                              % Bandwidth (KHz)
NoisePower    = 10^((N0+10*log10(B)+NF-30)/10);                     % dBm --> W
beta0         = (wavelength / (4*pi))^2;                            % average power gain at a reference distance of 1m

% IRS related parameters
d_irs_Y    = wavelength/2;                                          % length of the IRS element
d_irs_Z    = wavelength/2;                                          % width of the IRS element
H_i        = 25;                                                    % height of the IRS(m)
M          = [6 8 10 12 14 16 18 20];                               % M and N should be odd
N          = [6 8 10 12 14 16 18 20];
% M          = 6:30;                                                  % for testing power scaling law in test2.m
% N          = 6:30;
G_i        = [4 8];                                                 % antenna gain of IRS, 4->6dBi 8->9dBi
A          = 0.9;                                                   % amplitude of IRS element
alpha_y    = 0;                                                     % rotate angle along y-axis (vertical angle)
alpha_z    = 0;                                                     % rotate angle along z-axis (horizotal angle)

% TX related parameters
G_txdB     = 8;                                                     % antenna gain of TX,cited by 3GPP 38.901
G_tx       = 10^(G_txdB/10);
hpbw       = acosd( (1/2)^(1/(G_tx/2-1)) )*2;                       % half-power bandwidth,(cos hpbw/2)^(G_tx/2-1) = 1/2
P_t        = 10^((10-30)/10);                                       % 10dBm --> 0.01W
D          = zeros(1,length(M));                                    % TX-IRS distance(m)
for idx = 1:length(M)
    D(idx) = calDistanceFromTXtoIRS(M(idx)*N(idx), hpbw, d_irs_Y*d_irs_Z);
end

% UE related parameters
r             = 500/3;                                              % radius of sector
Hue           = [1.5 120];                                          % the height of UE(m) 1.5 for ground UE, 120 for aerial UE
numPlot       = length(Hue);                                        % num of UE layers for plotting

% Rician K-factor for aerial users
Kmin       = 0;                                                     % dB
Kmax       = 30;                                                    % dB
A1         = 10^(Kmin/10);                                          % environmental coefficient
A2         = 2/pi*log(10^(Kmax/10));

% Fixed BS pattern & 3D beamforming scheme
angleTiltV    = 100;                                                % vertical tilt angle
angleTiltH    = 0;                                                  % horizontal tilt angle
M_b           = [6 8 10 12 14 16 18 20];                            % num of vertical antennna elements
N_b           = [6 8 10 12 14 16 18 20];                            % num of horizontal antennna elements
% M_b          = 6:30;
% N_b          = 6:30;
dV            = wavelength/2;                                       % interval
dH            = wavelength/2;
centeridx     = floor(N_b./2).*M_b+floor(M_b./2)+1;                 % index of the center element

% Other parameters
numTrials       = 100;                                              % num of channel realizations
numfading       = 500;                                              % num of fading
numSectors      = 1;                                                % single cell
numBenchmark    = 2;                                                % num of the traditional BS scheme
numIRSscheme    = length(G_i);                                      % num of the IRS-aided BS scheme with different ERP
numScheme       = numBenchmark + numIRSscheme;

% angular spread from scattered waves
% cited by《A 3D geometry-based stochastic channel model for UAV-MIMO channels》
syms theta phi
k               = 0.5;                                              % spreading control parameter
theta_g         = pi/10;                                            % for MainIRSaidedBSinSingleCell/MultipleCells.m
theta_m         = pi/10; 
% theta_g         = [5*pi/12 pi/12];                                  % for RicianKfactorMCProve.m
% theta_m         = [pi/12 pi/12];                                    % for RicianKfactorMCProve.m
lower_cos       = pi/2 - (theta_g + theta_m);
upper_cos       = pi/2 - (theta_g - theta_m);
lower_von       = pi;
upper_von       = 2*pi;
phi_mu          = 3*pi/2;

% Fixed BS pattern & 3D beamforming scheme
theta_g_bs        = pi/10 * 180/pi;
theta_m_bs        = pi/10 * 180/pi;
lower_cos_bs      = 90 - (theta_g_bs + theta_m_bs);
upper_cos_bs      = 90 - (theta_g_bs - theta_m_bs);
theta_bs_ori      = (lower_cos_bs:1:upper_cos_bs)';

phi_bs_ori        = (-180:1:0)'; % range: -180~180
phi_mu_bs         = -90;
phi_bs            = reshape(repmat(phi_bs_ori', length(theta_bs_ori), 1), [], 1);
theta_bs          = repmat(90+theta_bs_ori, length(phi_bs_ori), 1); % ground

% Function
% cited by《3D trajectory optimization in Rician fading for UAV-enabled data harvesting》
funcRicianK = @(theta_K)A1.*exp(A2.*theta_K);
% cited by《A 3D geometry-based stochastic channel model for UAV-MIMO channels》
funcCosinePdf = @(theta,theta_g,theta_m)pi/(4*theta_m) * cos(pi/(2*theta_m)*(theta-theta_g));
funcVonMisesPdf = @(phi,phi_mu,k)exp(k*cos(phi-phi_mu))/(2*pi*besseli(0,k));

maxCosinePdf = pi./(4.*theta_m);
maxVonMisesPdf = exp(k)/(2*pi*besseli(0,k));

% Hypergeometric function
L = @(x)hypergeom(-1/2,1,x);         

% Figure Curve Color
color=[
    0    0.4470    0.7410      % deep blue 1
    0.8500    0.3250    0.0980 % orange 2
    0.9290    0.6940    0.1250 % deep yellow 3
    0.4940    0.1840    0.5560 % purple 4
    0.4660    0.6740    0.1880 % deep green 5
    0.3010    0.7450    0.9330 % sky blue 6
    0.8350    0.0780    0.1840 % wine red 7
    ];
