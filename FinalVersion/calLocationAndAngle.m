function [UE_2D_raw, UE_2D, sel_idx, dis_irs_tx_3D, b, theta_t, theta_tx] = calLocationAndAngle(numIRSIdx)
% CALLOCATIONANDANGLE calculate the location of TX, IRS, UE, departure
% angle of TX beam and incident angle of the IRS element
% Parameter explain:
% 1. numIRSIdx: index of the IRS element number
% 2. UE_2D_raw: 2D position of UE in a square region
% 3. UE_2D: 2D position of UE in a sector
% 4. sel_idx: index of select UEs in the sector
% 5. theta_t: incident angle of the IRS element
% 6. theta_tx: departure angle of the TX beam 
para_init;

%% coordinate
% 3D coordinate of IRS
% w0: the center of IRS
% w_I_bar: the coordinate of all IRS elements with respect to the IRS center
w0 = [0;0;25];
w_I_bar = [zeros(1, M(numIRSIdx)*N(numIRSIdx));
           reshape( repmat( ((1-M(numIRSIdx)/2:M(numIRSIdx)/2)-1/2).*d_irs_Y, N(numIRSIdx), 1 ), 1, [] );
           repmat( ((1-N(numIRSIdx)/2:N(numIRSIdx)/2)-1/2).*d_irs_Z, 1, M(numIRSIdx) )];
w_I = (w0 + w_I_bar)';

% 3D coordinate of TX
w_TX = [D(numIRSIdx) 0 H_i];

% 2D coordinate of UE
% create a square region
temp_ue_x = 0:3:2*r; % x-cord
temp_ue_y = -r*cosd(30):10:r*cosd(30); % y-cord
UE_x = reshape( repmat(temp_ue_x,length(temp_ue_y),1),[],1 );
UE_y = reshape( repmat(temp_ue_y',1,length(temp_ue_x)),[],1 );
UE_2D_raw = [UE_x UE_y];

% select UE in sector 1
k_upperLeft = temp_ue_y(end) / (0.5*r);
k_bottomLeft = temp_ue_y(1) / (0.5*r);
k_upperRight = (temp_ue_y(end)-0) / (1.5*r-2*r);
k_bottomRight = (temp_ue_y(1)-0) / (1.5*r-2*r);
b_upperRight = 0 - k_upperRight*2*r;
b_bottomRight = 0 - k_bottomRight*2*r;
sel_idx = find(UE_y - k_upperLeft*UE_x <=0 &...
               UE_y - k_bottomLeft*UE_x >=0 &...
               UE_y - k_upperRight*UE_x-b_upperRight <=0 &...
               UE_y - k_bottomRight*UE_x-b_bottomRight >=0 );
UE_2D = UE_2D_raw(sel_idx,:);

%% angle
% TX-IRS angle: theta_t and theta_tx
% b: the boresight of IRS element
% vecIRStoTx, vecTXtoIRScenter: the vector from IRS to TX, from TX to IRS center, respectively
% dis_irs_tx_3D, dis_tx_irs_center_3D: the 3D distance between TX and IRS, between TX and IRS center, respectively

b = [1 0 0];
vecIRStoTx = w_TX - w_I;
dis_irs_tx_3D = sqrt(sum(vecIRStoTx.^2,2));
theta_t = acosd( vecIRStoTx*b' ./ dis_irs_tx_3D );

vecTXtoIRScenter = w0'-w_TX;
dis_tx_irscenter_3D = sqrt(sum(vecTXtoIRScenter.^2,2));
theta_tx = abs(acosd( -vecIRStoTx*vecTXtoIRScenter' ./ (dis_irs_tx_3D.*dis_tx_irscenter_3D) ));

end
