function pulse = RFToPulse(RF, duration)
% RF is a complex valued vector.
% duration is in ms

numSteps = length(RF);

pulse.tp = duration; 
pulse.RFamp = abs(RF);
pulse.RFphase = phase(RF);
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);