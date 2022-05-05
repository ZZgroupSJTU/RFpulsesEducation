function pulse = PulseCreateGOIA(B1Max, duration, BW, sliceCenter, sliceThickness, B1ShapeNorm, gradShapeNorm, gradAxis)
% Create a GOIA frequency swept pulse, based on Tannus and Garwood, NMR in Biomed, 1997.
% Requires the normalized waveforms of the RF amplitude and gradient.
% Outputs a pulse structure.
%
% Input Variables
% Variable Name     Units     Description
% B1Max             kHz       Peak RF power
% duration          ms        Total pulse duration
% BW                kHz       Sweep bandwidth
% sliceCenter       mm        Self explanatory
% sliceThickness    mm        Self explanatory
% B1ShapeNorm       -         Normalized B1 amplitude waveform.
% gradShapeNorm     -         Normalized gradient waveform.
%
% Output Variables
% Variable Name     Units     Description
% pulse             -         Output pulse structure.
%                             pulse.tp       Pulse duration, in ms
%                             pulse.RFamp    Pulse amplitude, in kHz
%                             pulse.RFphase  Pulse phase, in radians
%                             pulse.Gx       x-gradient, in kHz/mm
%                             pulse.Gy       y-gradient, in kHz/mm
%                             pulse.Gz       z-gradient, in kHz/mm

if nargin<8, gradAxis = 'z'; end

numPoints      = numel(B1ShapeNorm);
dwellTime      = duration/numPoints; % ms

% Calculate the gradient amplitude, in kHz/mm
grad           = BW/sliceThickness;

% Calculate the RF's instantaneous frequency as a function of time (kHz*rad)
alpha          = cumsum(B1ShapeNorm.^2./gradShapeNorm*dwellTime); % ms
alpha0         = alpha(round(numPoints/2));
alpha          = alpha - alpha0; % Make sure wRF(t=Tp/2) = 0
wRFShapeNorm   = gradShapeNorm.*alpha;
wRFShapeNorm   = wRFShapeNorm./max(abs(wRFShapeNorm));
wRF            = wRFShapeNorm*2*pi*BW/2;

% Calculate the required additional phase to offset the pulse's center
% (A naive linear phase would distort the GOIA profile, not shift it, due
% to the time-dependent gradient)
wRFc           = 2*pi*gradShapeNorm*grad*sliceCenter; 
wRF            = wRF + wRFc;

% Calculate the RF pulse's phase, in radians
phi            = cumsum(wRF * dwellTime);

% Create the output pulse structure
pulse.tp       = duration;
pulse.RFamp    = B1ShapeNorm*B1Max;
pulse.RFphase  = phi;

pulse.Gx       = zeros(1,numPoints);
pulse.Gy       = zeros(1,numPoints);
pulse.Gz       = zeros(1,numPoints);
switch lower(gradAxis)
    case 'x'
        pulse.Gx       = grad*gradShapeNorm; % kHz/mm
    case 'y'
        pulse.Gy       = grad*gradShapeNorm; % kHz/mm
    case 'z'
        pulse.Gz       = grad*gradShapeNorm; % kHz/mm
    otherwise
        error('gradAxis = %s must be x, y, or z.', gradAxis);
end

