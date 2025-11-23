function AV=ElementPowerPatternVertical(angleV)
angle3dB=65;%%degree
Am=30;%%dB
AV=-min(12.*((angleV-90)./angle3dB).^2,Am);%%dB
