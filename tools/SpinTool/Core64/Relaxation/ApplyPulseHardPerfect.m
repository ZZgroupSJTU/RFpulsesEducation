function spins = ApplyPulseHardPerfect(spins, duration, flipAngle, pulsePhase)
% spins = ApplyPulseHard(duration, flipAngle, pulsePhase)
%
% Applies a hard pulse to a spin structure, ignoring all B1+ inhomogeneity
% Inputs 
%
% Name        Units         Description      
% Duration    microseconds  Pulse duration
% flipAngle   deg.          Desired flip angle
% pulsePhase  deg.          RF phase

pulse = PulseCreateHard(duration, flipAngle, pulsePhase);
spins = ApplyPulsePerfect(spins, pulse);

