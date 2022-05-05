function pulse = PulseCreateREBURP(duration)


An = [0.49 -1.02 1.11 -1.57 0.83 -0.42 0.26 -0.16 0.10 -0.07 0.04 -0.03 0.01 -0.02 0 -0.01]*1/duration;


numSteps = 256;
dt = duration/numSteps;
timeAxis = [0:dt:(numSteps-1)*dt];
w = 2*pi/duration;

RFEnvelope = zeros(1,numSteps);
for idx=1:numel(An)
    RFEnvelope = RFEnvelope + An(idx)*cos((idx-1)*w*timeAxis); 
end

pulse.tp = duration;
pulse.RFamp = RFEnvelope;
pulse.RFphase = zeros(1,numSteps);
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);
