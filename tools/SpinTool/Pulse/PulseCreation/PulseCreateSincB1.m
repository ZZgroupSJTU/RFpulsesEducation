function pulse = PulseCreateSincB1(numLobes, maxB1, numSteps, flipAngle)
% Creates a sinc small tip angle pulse with a given maximal B1, and
% calibrates its duration given the flip angle

% Convert flip angle to radians
flipAngle = flipAngle/180*pi;

timeAxis = linspace(-1/2, 1/2, numSteps);
dt = timeAxis(2)-timeAxis(1);
RFShape = maxB1.*sinc(timeAxis*numLobes);
duration = flipAngle/(2*pi*sum(RFShape*dt));


pulse.tp = duration;
pulse.RFamp = RFShape;
pulse.RFphase = zeros(1,numSteps);
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);

% % Calibrate for flip angle
% dt= duration/numSteps;
% RFIntegral = sum(pulse.RFamp*dt);  % int(B1*dt) = flipAngle
% pulse.RFamp = pulse.RFamp*flipAngle/RFIntegral/(2*pi);