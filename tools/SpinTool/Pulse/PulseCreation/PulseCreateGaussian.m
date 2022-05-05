function pulse = PulseCreateGaussian(centerFreq, widthFreq, SW, tiltAngle)
% The current function creates a pulse structure for an AM Gaussian pulse.
% The pulse's duration is calculated indirectly from its bandwidth, via
% Duration = 1/(Bandwidth*2.5).
%
% Input:
%
% Variable Name    Units         Description
% centerFreq       kHz           Center irradiation frequency
% widthFreq        kHz           FWHM of Gaussian excitation profile in 
%                                frequency domain.
% SW               kHz           Total spectral width of pulse (used to 
%                                determine the number of steps: 
%                                dwellTime=1/SW). Aliasing will occur 
%                                beyond SW limits.
% tiltAngle        degrees       Total tilt angle of pulse on resonance.

tiltAngle = tiltAngle/180*pi; % Convert to radians
dt = 1/SW; % ms

% Width in the time domain
widthTime = 1/widthFreq/2.5;
widthFactor = 4.2;
pulseDuration = widthFactor*widthTime;
numSteps = ceil(pulseDuration/dt);
tt = [0:dt:(numSteps-1)*dt];

% normalized
RF = exp(-(tt-pulseDuration/2).^2/widthTime^2);
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

