function pulseOut = PulseResample(pulseIn, N)



tAxisIn = linspace(0, pulseIn.tp, numel(pulseIn.RFamp));
tAxisOut = linspace(0, pulseIn.tp, N);

pulseOut.tp      = pulseIn.tp;
pulseOut.RFamp   = interp1(tAxisIn, pulseIn.RFamp, tAxisOut);
pulseOut.RFphase = interp1(tAxisIn, pulseIn.RFphase, tAxisOut);
pulseOut.Gx      = interp1(tAxisIn, pulseIn.Gx, tAxisOut);
pulseOut.Gy      = interp1(tAxisIn, pulseIn.Gy, tAxisOut);
pulseOut.Gz      = interp1(tAxisIn, pulseIn.Gz, tAxisOut);

