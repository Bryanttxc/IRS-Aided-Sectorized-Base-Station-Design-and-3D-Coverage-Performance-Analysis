function [GarraydB]=ArrayPowerPattern(angleV,angleH,angleTiltV,angleTiltH,M,N,dV,dH,wavelength)

GelementdB=ElementPowerPatternOverall(angleV,angleH); % dBi
Gelement=10.^(GelementdB./10); % 转换公式 dBi --> W
Felement=sqrt(Gelement); % signal

k=2*pi/wavelength;% wave number
kVector=k.*[sin(angleV./180*pi).*cos(angleH./180*pi),sin(angleV./180*pi).*sin(angleH./180*pi),cos(angleV./180*pi)];%wave vector 极坐标转换为[x,y,z]直角坐标
rMatrix=zeros(M*N,3);%%location of the antenna element (x,y,z) pile up (into a long matrix), reference origin is the lower left element
cnt = 1;
idx = -N/2+1:N/2;
for n=1:N
     rMatrix(1+(cnt-1)*M:cnt*M,3) = [-M/2:(M/2-1)]'*dV;
     if N == 1
        rMatrix(1+(cnt-1)*M:cnt*M,2) = (n-1)*dH;
     elseif mod(N,2)==0
        rMatrix(1+(cnt-1)*M:cnt*M,2) = (idx(n)-1)*dH;        
     end
     cnt = cnt + 1;
end
SteeringVector=exp(-1i.*(rMatrix*kVector')); %导向矢量 位置*入射角
w=1./sqrt(M).*exp(-1i.*k.*([-M/2:(M/2-1)]').*dV.*cos(angleTiltV/180*pi));%%vertical weight vector, M-by-1, with unit power, to form angleTiltV
v=1./sqrt(N).*exp(-1i.*k.*([0:(N-1)]).*dH.*sin(angleTiltH/180*pi));%%horizontal weight vector, 1-by-N, to form angleTiltH
Weight2D=kron(v,w);%%M-by-N weight matrix
WeightFlatten=reshape(Weight2D,[M*N,1]); %%M*N-by-1 weight vector, taken column-wise from Weight2D
ArrayFactor=WeightFlatten'*SteeringVector;% 加权
Farray=Felement.*ArrayFactor';
Garray=(abs(Farray)).^2;
GarraydB=10.*log10(Garray); % ULA天线增益 W-->dB
