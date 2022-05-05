function spinsOut = PurgeMoment3D(spinsIn, kx, ky, kz, tp)
% SYNTAX: spinsOut = PurgeMoment(spinsIn, kx, ky, kz, tp)
%
% Applies a constant gradient (mT/m) to the spins for a duration tp (ms).
% This takes into account relaxation. The moments kx, ky, kz are given
% in meters^(-1).
%
% k/tp ~ m^(-1)/ms = kHz/m 
% 0.001*k/tp = kHz/mm

if tp==0
    spinsOut = spinsIn;
else
    Gx = 0.001*kx/tp; % kHz/mm
    Gy = 0.001*ky/tp; % kHz/mm
    Gz = 0.001*kz/tp; % kHz/mm
    pulse.tp      = tp;
    pulse.RFamp   = 0; 
    pulse.RFphase = 0;
    pulse.Gx      = Gx;
    pulse.Gy      = Gy;
    pulse.Gz      = Gz;
    spinsOut = ApplyPulse3D(spinsIn,pulse);
end