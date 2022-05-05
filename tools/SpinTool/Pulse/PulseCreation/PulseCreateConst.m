function pulse = PulseCreateConst(duration,numSteps,pulsePower,pulsePhase)
% Creates a constant RF pulse structure.
% duration     ms
% numSteps   
% pulsePower   kHz
% pulsePhase   radians

pulse.tp      = duration;
pulse.RFamp   = ones(1,numSteps)*pulsePower;
pulse.RFphase = ones(1,numSteps)*pulsePhase;
pulse.Gx      = zeros(1,numSteps);
pulse.Gy      = zeros(1,numSteps);
pulse.Gz      = zeros(1,numSteps);