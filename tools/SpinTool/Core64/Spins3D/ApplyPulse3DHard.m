function spins = ApplyPulse3DHard(spins, duration, flipAngle, pulsePhase)
% spins = ApplyPulse3DHard(duration, flipAngle, pulsePhase)
%
% Applies a hard pulse to a spin structure.
% Inputs 
%
% Name        Units         Description      
% Duration    microseconds  Pulse duration
% flipAngle   deg.          Desired flip angle
% pulsePhase  deg.          RF phase

pulse = PulseCreateHard(duration, flipAngle, 0);
pulsePhase = pulsePhase/180*pi;
spins = ApplyPulse3DCycle(spins, pulse, pulsePhase);