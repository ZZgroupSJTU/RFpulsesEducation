function pulse = PulseAddPhase(pulse, ph)
% SYNTAX: function pulse = PulseAddPhase(pulse, ph)
%
% Add a phase ph to the pulse.

pulse.RFphase = pulse.RFphase + ph;