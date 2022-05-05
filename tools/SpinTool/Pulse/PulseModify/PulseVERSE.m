function pulseVERSE = PulseVERSE(pulse, modFun, gradAxis)
% SYNTAX:
%
%     pulse = PulseVERSE(pulse, modFun, [gradAxis])
%
% Returns the VERSE-ified version of the input pulse given the gradient
% modulation function modFun.  
%
% Example: 
%     % Create a sinc pulse
%     N = 512;
%     pulse = PulseCreateSinc(4, 12.8, N, 90);
%     G0 = 0.2;
%     pulse.Gz = G0*ones(1,N);
%     % Create the modulation function (don't worry about total area)
%     modFun = (1 - 0.8*sin([0:N-1]/(N-1)*pi));
%     % VERSE-ify it.
%     pulseVERSE = PulseVERSE(pulse, modFun);
%
% Based on Conolley et. al., MRM 18, 28-38 (1991).
%
% Notes:
% [1] The modulation function should integrate to T, the total duration of 
%     the pulse. If it doesn't, its amplitude will be adjusted until its 
%     area matches T.
% [2] The duration of the VERSE-ified pulse will be equal to that of the 
%     input pulse.
% [3] To further reduce the peak amplitude of the pulse and gradient at
%     the expense of gradient duration, perform the transformation:
%     pulse.RFamp -> pulse.RFamp/a
%     pulse.Gz    -> pulse.Gz/a
%     pulse.tp    -> pulse.tp*a
% [4] In general, for a sinc like pulse you'd like modFun to dip at the
%     center and be high at the edges, which will lead to the RF amplitude
%     becoming higher at the edges (where it's low to begin with) and
%     smaller at the center (where its peak is high).
% [5] Don't forget: this produces bad off resonance behavior! Unlike for
%     regular constant gradient pulses, in which a chemical shift produces
%     a proportional spatial shift, here the entire profile gets distorted.

N = numel(pulse.RFamp);
dt = pulse.tp/N;
timeAxisOrig = [0:dt:pulse.tp-dt];
timeAxis = [0:dt:pulse.tp];

if nargin<3
    gradAxis = 'z';
end

RFamp = pulse.RFamp;
RFphase = pulse.RFphase;

switch lower(gradAxis)
    case 'x'
        gradShape = pulse.Gx;
    case 'y'
        gradShape = pulse.Gy;
    case 'z'
        gradShape = pulse.Gz;
    otherwise
        error('gradAxis has to be x, y or z.');
end

% To ensure the gridding takes place successfully, we need to extrapolate
% the RF and gradient waveform values all the way to pulse.tp (they only
% reach until pulse.tp-dt)
RFamp(end+1) = RFamp(end);
RFphase(end+1) = RFphase(end);
gradShape(end+1) = gradShape(end);

timeAxisModFun = linspace(0, pulse.tp, numel(modFun));

% Interpolate the modulation and time axis to high resolution
OSFactor = 20;
dtNew = pulse.tp/(N*OSFactor);
timeAxisNew = [0:dtNew:pulse.tp];
modFunNew = interp1(timeAxisModFun, modFun, timeAxisNew);
modFunNew = modFunNew.*(pulse.tp/sum(modFunNew.*dtNew)); % Make sure it integrates to pulse.tp
RFAmpNew = interp1(timeAxis, RFamp, timeAxisNew);
RFPhaseNew = interp1(timeAxis, RFphase, timeAxisNew);
gradNew = interp1(timeAxis, gradShape, timeAxisNew);

% Calculate the non-uniform time axis for VERSE
timeAxisVERSE = [0, cumsum(modFunNew(1:end-1)*dtNew)];

% Resample RF and gradient waveforms on non-uniform grid
RFAmpVERSE = interp1(timeAxisNew, RFAmpNew, timeAxisVERSE);
RFPhaseVERSE = interp1(timeAxisNew, RFPhaseNew, timeAxisVERSE);
gradVERSE = interp1(timeAxisNew, gradNew, timeAxisVERSE);

% Downsample to the same number of points as the original pulse
pulseVERSE.tp = pulse.tp;
pulseVERSE.RFamp = interp1(timeAxisNew, RFAmpVERSE.*modFunNew, timeAxisOrig);
pulseVERSE.RFphase = interp1(timeAxisNew, RFPhaseVERSE, timeAxisOrig);
gradVERSEInterp = interp1(timeAxisNew, gradVERSE.*modFunNew, timeAxisOrig);
pulseVERSE.Gx = zeros(1, N);
pulseVERSE.Gy = zeros(1, N);
pulseVERSE.Gz = zeros(1, N);
switch lower(gradAxis)
    case 'x'
        pulseVERSE.Gx = gradVERSEInterp;
    case 'y'
        pulseVERSE.Gy = gradVERSEInterp;
    case 'z'
        pulseVERSE.Gz = gradVERSEInterp;
end


