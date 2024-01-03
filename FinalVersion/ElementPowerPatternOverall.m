function GelementdB=ElementPowerPatternOverall(angleV,angleH)
Gemax=8;%%dBi 功率增益单位 10dBi=10lg(10)
Am=30;%%dB
Aelement=-min(-(ElementPowerPatternVertical(angleV)+ElementPowerPatternHorizontal(angleH)),Am);%%dB
GelementdB=Gemax+Aelement;%dBi