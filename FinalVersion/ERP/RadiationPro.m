function [F] = RadiationPro(theta,label,G)
% Previous version of getRadiationValue.m

[x,y] = size(theta);
F = zeros(x,y);

if strcmp(label,'tx') == 1 || strcmp(label,'rx') == 1
    for p = 1:x
        for q = 1:y
            if theta(p,q) <= 90 && theta(p,q) >= 0
                F(p,q) = cosd(theta(p,q)).^(G/2-1);
            elseif theta(p,q) > 90 && theta(p,q) <= 180
                F(p,q) = 0;
            else
                sprintf("the value of theta is error!")
            end
        end
    end
elseif strcmp(label,'t') == 1 || strcmp(label,'r') == 1
    for p = 1:x
        for q = 1:y
            if theta(p,q) <= 90 && theta(p,q) >= 0
                F(p,q) = cosd(theta(p,q))^(G/2-1);           
            elseif theta(p,q) > 90 && theta(p,q) <= 180
                F(p,q) = 0;
            else
                sprintf("the value of theta is error!")
            end
        end
    end
end


end