function pulse = PulseCreateGaussianSiemens(centerFreq, pulseDuration, numSteps, tiltAngle)
% The current function creates a pulse structure for an AM Gaussian pulse.
% The pulse's duration is calculated indirectly from its bandwidth, to
% match the creation procedure on a Siemens machine.
%
% Input:
%
% Variable Name    Units         Description
% centerFreq       kHz           Center irradiation frequency
% pulseDuration    ms            Pulse duration
% numSteps        -              Number of pulse points
% tiltAngle        degrees       Total tilt angle of pulse on resonance.


tiltAngle = tiltAngle/180*pi; % Convert to radians
dt = pulseDuration/numSteps; % ms
tt = [0:dt:(numSteps-1)*dt];
widthTime = pulseDuration*0.35;

% normalized
RF = exp(-(tt-pulseDuration/2).^2/widthTime^2);
RF = RF-RF(1);

RFAngle = (2*pi*sum(RF*dt));

% RFTiltAngle = (2*pi*sum(RF*dt)); % On resonance
maxB1 = tiltAngle/RFAngle;

% Create time vector
pulse.tp = pulseDuration;
pulse.RFamp = RF*maxB1;
pulse.RFphase = centerFreq*2*pi*tt;
pulse.Gx = 0.*RF;
pulse.Gy = 0.*RF;
pulse.Gz = 0.*RF;

