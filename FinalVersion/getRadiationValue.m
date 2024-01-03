function [F] = getRadiationValue(theta, G)
% GETRADIATIONVALUE applies for antennas whose radiation pattern expression is cos^(a)
% Parameter explain:
% Input
% 1. theta: denote the elevation angle respectively in the local coordinate 
%           system of the antenna, 0-180°;
% 2. G: maximum radiation gain corresponding to theta=0°.
% Output
% 3. F = cos^(G/2-1)(theta);
%   denotes the radiation pattern value of antennas.
%%%%%%%%%%%%%%%%%%%%%%%%
[x,y] = size(theta); % dimension of theta
F = zeros(x,y);

% index
idx0to90 = theta >= 0 & theta <= 90;
idx90to180 = theta > 90 & theta <= 180;

F(idx0to90 == 1) = cosd(theta(idx0to90 == 1)).^(G/2-1);
F(idx90to180 == 1) = 0;
F(idx0to90 == 0 & idx90to180 == 0) = -1;

% check if F has error values
errorIdx = find(F == -1);
if(~isnan(errorIdx))
    error("the value of theta is error!")
end

end
