function PrLOS = LOSprobability(Hue,d2D) %%only for a single UE, not for a batch of UEs
Hthreshold=100;
if Hue<=22.5 %Uma
    if d2D<=18
        PrLOS=1;
    else
        Cprim=0;
        if Hue>13
            Cprim=((Hue-13)/10)^1.5;
        end
        PrLOS=( 18/d2D+(1-18/d2D)*exp(-d2D/63) )*( 1+Cprim*5/4*(d2D/100)^3*exp(-d2D/150) );
    end
elseif Hue<=Hthreshold
    d1=max(18, 460*log10(Hue)-700);  % 36.4m
    p1=4300*log10(Hue)-3800;
    if d2D<=d1
        PrLOS=1;
    else
        PrLOS=d1/d2D+(1-d1/d2D)*exp(-d2D/p1);
    end
else
    PrLOS=1;
end
