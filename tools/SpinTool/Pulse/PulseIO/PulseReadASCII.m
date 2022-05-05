function pulse = PulseReadASCII(filename, duration, peakB1)
% SYNTAX:
%
%    pulse = PulseReadASCII(filename, duration, peakB1)
%
% 


pulseData = load(filename);
numSteps = size(pulseData,1);
pulse.tp = duration;
pulse.RFamp = pulseData(:,1)'*42.57; % kHz
pulse.RFphase = pulseData(:,2)'/360*2*pi; 
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);
