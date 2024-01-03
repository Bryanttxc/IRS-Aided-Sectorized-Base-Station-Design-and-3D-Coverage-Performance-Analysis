function PL_LOS=LOSpathloss(Hue,d2D,d3D,Hbs,fc,c)%%only for a single UE, not for a batch of UEs. in dB
HE=1;%m
HuePrim=Hue-HE;
HBprim=Hbs-HE;
dBP=4*HBprim*HuePrim*fc/c;
Hthreshold=300;
fcGHz=fc*10^(-9);
if Hue<=22.5
    if d2D<=dBP%%%minimal 10 meters in 3GPP TR38.901, not implemented here
        PL_LOS=28+22*log10(d3D)+20*log10(fcGHz);
    else%%%maximal 5km in 3GPP TR38.901, not implemented here
        PL_LOS=28+40*log10(d3D)+20*log10(fcGHz)-9*log10( dBP^2+(Hbs-Hue)^2 );        
    end
elseif Hue<=Hthreshold
    PL_LOS=28+22*log10(d3D)+20*log10(fcGHz);%%%maximal d2D=4km in 3GPP RAN1 #90, R1-1714856, not implemented here
else
    msg = 'Current 3GPP model only defines Aerial UE up to 300 meters in altitude.';
    error(msg)
end


