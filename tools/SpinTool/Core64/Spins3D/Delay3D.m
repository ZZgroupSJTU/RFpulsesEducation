function spinsOut = Delay3D(spinsIn, tp)
% SYNTAX: spinsOut = Delay3D(spinsIn, tp)
%
% Applies a delay, taking relaxation into account.

pulse.tp      = tp;
pulse.RFamp   = 0; 
pulse.RFphase = 0;
pulse.Gx      = 0;
pulse.Gy      = 0;
pulse.Gz      = 0;

% temp1 & temp2 aren't used
spinsOut = ApplyPulse3D(spinsIn,pulse);