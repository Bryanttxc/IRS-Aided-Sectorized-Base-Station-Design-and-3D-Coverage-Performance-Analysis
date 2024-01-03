function PL_NLOS=NLOSpathloss(Hue,d2D,d3D,Hbs,fc,c)%%only for a single UE, not for a batch of UEs. in dB
HE=1;%m
HuePrim=Hue-HE;
HBprim=Hbs-HE;
dBP=4*HBprim*HuePrim*fc/c;
Hthreshold=300;
fcGHz=fc*10^(-9);
if Hue<=22.5
    PL_LOS=LOSpathloss(Hue,d2D,d3D,Hbs,fc,c);
    PLprim=13.54+39.08*log10(d3D)+20*log10(fcGHz)-0.6*(Hue-1.5);
    PL_NLOS=max(PL_LOS,PLprim);
elseif Hue<=Hthreshold
    PL_NLOS=22.5+( 46-7*log10(Hue) )*log10(d3D) + 20*log10(40*pi*fcGHz/3);%%%maximal d2D=4km in 3GPP RAN1 #90, R1-1714856, not implemented here
else
    msg = 'Current 3GPP model only defines Aerial UE up to 300 meters in altitude.';
    error(msg)
end


