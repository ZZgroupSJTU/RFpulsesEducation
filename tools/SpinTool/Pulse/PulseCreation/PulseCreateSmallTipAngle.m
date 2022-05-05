function pulse = PulseCreateSmallTipAngle(freqAxis, freqResponse, T, B1)
% Generates a small-tip-angle approximate RF pulse using an inverse Fourier
% transform.
%
%   pulse = PulseCreateSmallTipAngle(freqAxis, freqResponse, maxB1)
%
% Input Variable
% Name           Units    Description
% freqAxis       kHz      Frequency axis of response
% freqResponse   -        Desired response (both mag. & phase)
% maxB1          kHz      Available peak power of B1.
%
% Output Variable
% pulse          -        Output pulse structure
%
% Theory: for small tip angles, 
%
%    Mxy(T) = i*M0*gamma*exp(-i*w*T/2)*Integrate[B1(t)*exp(i*w*t)*dt]
%           = i*M0*gamma*exp(-i*w*T/2)*FFT(B1(t))
%
% where B1(t) is the complex RF waveform, centered about t=0, and
% integration is carried out between -T/2 and T/2. By inverting this,
% we can approximate (taking M0=1):
%
%    B1(T) = IFFT[-i/gamma*exp(i*w*T/2)*Mxy(T)]

% ------------------------------------------------------------------------
% Step I: Interpolate in the frequency domain
%
% Interpolation serves to (1.) create a uniform sampling axis, which can
% be FFT-ed, and (2.) make dv smaller, thus reducing aliasing artifacts
% in the time domain (which we will be truncating).
% ------------------------------------------------------------------------

interpolationFactor = 10;

% Find the minimal dv in the provided frequency axis.
dv = min(abs(diff(freqAxis)));

% Make it smaller by an order of magnitude
dv = dv/interpolationFactor;

% Compute spectral width of given desired profile
SW = freqAxis(end) - freqAxis(1);

% Created interpolated frequency axis
NI = ceil(SW/dv);
freqAxisI = linspace(freqAxis(1), freqAxis(end), NI);

% Interpolate
freqResponseI = interp1(freqAxis, freqResponse, freqAxisI);

% ------------------------------------------------------------------------
% Step II: Extrapolate by zero-filling in the frequency domain
%
% Extrapolation (in the frequency domain) serves to interpolate in the 
% time domain, making the RF waveform smoother (and "more correct") before 
% resampling.
% ------------------------------------------------------------------------

extrapolationFactor = 3;
zerosVec = zeros(1, NI*extrapolationFactor);
freqResponseIE = [zerosVec, freqResponseI, zerosVec];
SWIE = SW*(extrapolationFactor*2+1);
NIE = numel(freqResponseIE);

% ------------------------------------------------------------------------
% Step III: FFT
% ------------------------------------------------------------------------

pulseShapeIE = fftshift(fft(fftshift(freqResponseIE)));
timeAxisIE = linspace(-SWIE/2, SWIE/2, NIE);

% ------------------------------------------------------------------------
% Step IV: Resample time-domain RF shape
% ------------------------------------------------------------------------

% Define the resampled time axis
dwellTime = 1/SW;
numSteps = ceil(T/dwellTime);
timeAxis = linspace(-T/2, T/2, numSteps);

% Resample
pulseShape = interp1(timeAxisIE, pulseShapeIE, timeAxis);

% ------------------------------------------------------------------------
% Step V: Create the desired pulse shape
% ------------------------------------------------------------------------

pulse.tp = T;
pulse.RFamp = abs(pulseShape)*B1;
pulse.RFphase = phase(pulseShape);
pulse.Gx = zeros(1, numSteps);
pulse.Gy = zeros(1, numSteps);
pulse.Gz = zeros(1, numSteps);