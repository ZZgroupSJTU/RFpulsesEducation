function pulseout = PulseToPositive(pulsein)
% SYNTAX:
%
% The current function makes sure pulse.RFamp is positive.
% Negative values of pulse.RFamp correspond to a phase flip of 180 degrees,
% so, e.g. RFamp = -1 and RFphase = pi describe a pulse step with an amplitude
% of 1 and a phase of 0. This function "corrects" the pulse such that, e.g.
% (in our example) RFamp = 1 and RFphase = 0.

N = length(pulsein.RFamp);
pulseout = pulsein;

% Recall sgn = +1 for positive, -1 for negative, and 0 for 0.
sgn = sign(pulsein.RFamp);
% make sure sgn(i) = 1 when it is negative, 0 otherwise
sgn = (sgn<0);

% For negative RFamp, invert RFamp and add a pi phase to the pulse.
pulseout.RFamp = abs(pulseout.RFamp);
pulseout.RFphase = pulseout.RFphase + sgn*pi;

% Make sure the phase of the resulting pulse is between [0, 2pi]
pulseout.RFphase = mod(pulseout.RFphase, 2*pi);