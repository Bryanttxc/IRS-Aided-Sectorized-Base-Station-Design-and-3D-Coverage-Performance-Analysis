% 水平方向上的能量模式
function AH=ElementPowerPatternHorizontal(angleH)
angle3dB=65;%%degree
Am=30;%%dB
AH=-min(12.*(angleH./angle3dB).^2,Am);%%dB 
%Explanation:the above formula is derived from the report of 3GPP
