function pulse = PulseConcat(pulse1,pulse2);
% SYNTAX: function pulse = PulseConcat(pulse1,pulse2);
%
% Concatenates the two pulses, pulse1 & pulse2.
% Assumption: dwell time is the same for both.

pulse.tp      = pulse1.tp + pulse2.tp;
pulse.RFamp   = [pulse1.RFamp,pulse2.RFamp];
pulse.RFphase = [pulse1.RFphase,pulse2.RFphase];
pulse.Gx      = [pulse1.Gx,pulse2.Gx];
pulse.Gy      = [pulse1.Gy,pulse2.Gy];
pulse.Gz      = [pulse1.Gz,pulse2.Gz];
