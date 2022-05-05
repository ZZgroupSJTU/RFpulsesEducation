function rfstruct=calc_SLRrf(rfstruct,theta,phi,refpw,reftpwr)

% n = rfstruct.pts;
% polyn = length(theta);
% theta = theta(polyn/2-n/2+1:polyn/2+n/2);
% phi = phi(polyn/2-n/2+1:polyn/2+n/2);

rfstruct.res=rfstruct.duration/rfstruct.pts;
B1amp=theta / (1e3*2*pi*rfstruct.res);              %[kHz]

%--------------calc_SLRrf power-------------------------%
refB1amp=1e-3/(4*refpw);
maxB1amp=max(B1amp);
% rfstruct.B1max=maxB1amp;
rfpwr= reftpwr+20*log10(maxB1amp/refB1amp);
rfpwr1=4095*10^(-(ceil(rfpwr)-rfpwr)/20);
rfstruct.powerCoarse=ceil(rfpwr);
rfstruct.powerFine=round(rfpwr1);

%---------------updata rfstruct-------------------------%

rfstruct.phase = mod(phi*180/pi,360);
rfstruct.amp=B1amp;

t=(0.5*rfstruct.res)+linspace(0,rfstruct.duration-rfstruct.res,rfstruct.pts);
rfstruct.head.integral=trapz(t,B1amp);

rfstruct.head.bandwidth     =rfstruct.bandwidth;

rfstruct.head.inversionBw   =0;  
rfstruct.head.modulation    ='phase';
rfstruct.head.rfFraction    =0;
rfstruct.head.type          ='SLR';
rfstruct.head.version       =1.0;