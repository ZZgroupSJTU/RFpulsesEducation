function pulse = PulseCreateGaussianTime(centerFreq, pulseDuration, widthTime, SW, tiltAngle)
% The current function creates a pulse structure for an AM Gaussian pulse.
% Unlike the function PulseCreateGaussian, this function accepts the
% duration of the Gaussian, not its frequency-domain bandwidth.
%
% Input:
%
% Variable Name    Units         Description
% centerFreq       kHz           Center irradiation frequency
% pulseDuration    ms            Duration of pulse.
% widthTime        ms            FWHM of pulse in time domain
% SW               kHz           Total spectral width of pulse (used to 
%                                determine the number of steps: 
%                                dwellTime=1/SW). Aliasing will occur 
%                                beyond SW limits.
% tiltAngle        degrees       Total tilt angle of pulse on resonance.

tiltAngle = tiltAngle/360*2*pi; % Convert to radians
dt = 1/SW; % ms
numSteps = ceil(pulseDuration/dt);
tt = [0:dt:(numSteps-1)*dt];

% normalized
RF = exp(-(tt-pulseDuration/2).^2/widthTime^2);
RF = RF - RF(1);

maxB1 = tiltAngle/sum(2*pi*RF*dt);

% Create time vector
pulse.tp = pulseDuration;
pulse.RFamp = maxB1*RF;
pulse.RFphase = centerFreq*2*pi*tt;
pulse.Gx = 0.*RF;
pulse.Gy = 0.*RF;
pulse.Gz = 0.*RF;