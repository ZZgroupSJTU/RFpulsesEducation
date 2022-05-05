function isPulseAM = IsPulseAM(pulse, errorThreshold)
% IsPulseAM  Returns true if pulse is amplitude-modulated

if nargin<2, errorThreshold = 1e-4; end
pulsePhase = pulse.RFphase;
pulsePhase = mod(pulsePhase, 2*pi);
numStepsPulse = numel(pulsePhase);
numStepsPIPhase = numel(find((abs(pulsePhase-pi)<errorThreshold) + (abs(pulsePhase)<errorThreshold)));
isPulseAM = (numStepsPulse == numStepsPIPhase);

