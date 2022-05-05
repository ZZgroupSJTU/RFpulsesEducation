function pulse = PulseCreateHyperGaussian(centerFreq, widthFreq, SW, tiltAngle)
% The current function creates a pulse structure for an AM Hyper-Gaussian 
% pulse, i.e. one with a profile exp(-x^4) (instead of exp(-x^2) for a
% regular Gaussian).
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
b = 0.5;
widthTime = 1/widthFreq/2.5*b;
widthFactor = 4.2/b;
pulseDuration = widthFactor*widthTime;
numSteps = ceil(pulseDuration/dt);
tt = [0:dt:(numSteps-1)*dt];
tt = tt - pulseDuration/2; 

% normalized
RF = exp(-tt.^2/widthTime^2).*sinc(tt/2);

% Trim edges
a = 0.3;
RF = RF(ceil(a*numSteps):floor(numSteps*(1-a))).*(1-2*a);
tt = tt(ceil(a*numSteps):floor(numSteps*(1-a)));

% Create time vector
pulse.tp = pulseDuration*2;
pulse.RFamp = RF;
pulse.RFphase = centerFreq*2*pi*tt;
pulse.Gx = 0.*RF;
pulse.Gy = 0.*RF;
pulse.Gz = 0.*RF;

% Calibrate power
numSteps = numel(pulse.RFamp);
dt = pulse.tp/numSteps;
RFAngle = (2*pi*sum(pulse.RFamp*dt))
maxB1 = tiltAngle/RFAngle
pulse.RFamp = pulse.RFamp*maxB1
