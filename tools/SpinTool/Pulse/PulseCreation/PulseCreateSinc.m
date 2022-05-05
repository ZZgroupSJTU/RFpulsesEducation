function pulse = PulseCreateSinc(numLobes, duration, numSteps, flipAngle, pulsePhase)
% SYNTAX:
%
%    pulse = PulseCreateSinc(numLobes, duration, numSteps, flipAngle, [pulsePhase in deg.])
%
% Creates a sinc pulse.

if nargin<5
    pulsePhase = 0;
end

% Convert flip angle to radians
flipAngle = flipAngle/180*pi;

timeAxis = linspace(-numLobes, numLobes, numSteps);

pulse.tp = duration;
pulse.RFamp = sinc(timeAxis);
pulse.RFphase = ones(1,numSteps)*pulsePhase/180*pi;
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);

% Calibrate for flip angle
dt= duration/numSteps;
RFIntegral = sum(pulse.RFamp*dt);  % int(B1*dt) = flipAngle
pulse.RFamp = pulse.RFamp*flipAngle/RFIntegral/(2*pi);