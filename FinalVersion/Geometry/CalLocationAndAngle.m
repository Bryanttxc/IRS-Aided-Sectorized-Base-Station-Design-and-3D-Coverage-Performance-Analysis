function [w_UE_2D_squ, w_UE_2D, sel_idx, dis_IRS2TX_3D, b, theta_t, theta_tx] = CalLocationAndAngle(num_IRS_idx)
% CALLOCATIONANDANGLE calculate the position of TX, IRS, UE, departure
% angle of TX beam and incident angle of the IRS element
% Parameter explain:
% 1. num_IRS_idx: index of the IRS element number
% 2. w_UE_2D_squ: 2D position of UE in a square region
% 3. w_UE_2D: 2D position of UE in a sector
% 4. sel_idx: index of select UEs in the sector
% 5. theta_t: incident angle of the IRS element
% 6. theta_tx: departure angle of the TX beam 

para_init;

% handle func to gen coordinates of odd/even elements
hdfunc_cord_gen = @(num_elem)(1-floor(num_elem/2):floor(num_elem/2))'-(1/2*(1-mod(num_elem,2)));
hdfunc_angle_gen = @(vec_a, vec_b, abs_a, abs_b)acosd( vec_a*vec_b' ./ (abs_a*abs_b));

%% 1) generate coordinates
% 3D coordinate of IRS
% w0: the center of IRS
% w_I_bar: the coordinates of all IRS elements with respect to the IRS center
w0 = [0 0 H_i];
tot_IRS_elem = M(num_IRS_idx)*N(num_IRS_idx);
x_I_bar = zeros(tot_IRS_elem, 1);
y_I_bar = reshape( repmat( hdfunc_cord_gen(M(num_IRS_idx)).*d_IRS_Y, 1, N(num_IRS_idx))', [], 1);
z_I_bar = repmat( hdfunc_cord_gen(N(num_IRS_idx)).*d_IRS_Z, M(num_IRS_idx), 1);

w_I_bar = [x_I_bar y_I_bar z_I_bar];
w_I = w0 + w_I_bar;

% 3D coordinates of TX, same height as IRS
w_TX = [D(num_IRS_idx) 0 H_i];

% 2D coordinates of UE
% create a square region
% ----------------
% | /          \ |
% |              |
% |   sector 1   |
% |              |
% | \          / |
% ----------------
x_ue = 0 : sample_interv_x : 2*r;
y_ue = -r*cosd(30) : sample_interv_y : r*cosd(30);
x_UE = reshape( repmat(x_ue,length(y_ue),1),[],1 );
y_UE = reshape( repmat(y_ue',1,length(x_ue)),[],1 );

w_UE_2D_squ = [x_UE y_UE];

% select UE in sector 1
%   ------------
%   /          \
% |              |
% |   sector 1   |
% |              |
%   \          /
%   ------------
k_upper_left = y_ue(end) / (0.5*r);
k_bottom_left = y_ue(1) / (0.5*r);
k_upper_right = (y_ue(end)-0) / (1.5*r-2*r);
k_bottom_right = (y_ue(1)-0) / (1.5*r-2*r);
b_upper_right = 0 - k_upper_right*2*r;
b_bottom_right = 0 - k_bottom_right*2*r;
sel_idx = find(y_UE - k_upper_left*x_UE <=0 &...
               y_UE - k_bottom_left*x_UE >=0 &...
               y_UE - k_upper_right*x_UE-b_upper_right <=0 &...
               y_UE - k_bottom_right*x_UE-b_bottom_right >=0 );

w_UE_2D = w_UE_2D_squ(sel_idx,:);

%% 2) generate angles
% TX-IRS angle: theta_t and theta_tx
% b: the boresight of IRS element
b = [1 0 0];
vec_IRS2TX = w_TX - w_I;
dis_IRS2TX_3D = sqrt(sum(vec_IRS2TX.^2,2));
theta_t = hdfunc_angle_gen(vec_IRS2TX, b, dis_IRS2TX_3D, 1);

vec_TX2IRS = -vec_IRS2TX;
dis_TX2IRS_3D = dis_IRS2TX_3D;
vec_TX2IRScenter = w0 - w_TX;
dis_TX2IRScenter_3D = sqrt(sum(vec_TX2IRScenter.^2,2));
theta_tx = abs(hdfunc_angle_gen(vec_TX2IRS, vec_TX2IRScenter, dis_TX2IRS_3D, dis_TX2IRScenter_3D));

end
