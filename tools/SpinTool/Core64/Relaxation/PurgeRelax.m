function spinsOut = PurgeRelax(spinsIn, Gx, Gy, Gz, tp)
% SYNTAX: spinsOut = PurgeRelax(spinsIn, Gx, Gy, Gz, tp)
%
% Applies a constant gradient (kHz/mm) to the spins for a duration tp (ms).
% This takes into account relaxation

pulse.tp      = tp;
pulse.RFamp   = [0]; 
pulse.RFphase = [0];
pulse.Gx      = [Gx];
pulse.Gy      = [Gy];
pulse.Gz      = [Gz];

% temp1 & temp2 aren't used
spinsOut = ApplyPulseRelax(spinsIn,pulse);