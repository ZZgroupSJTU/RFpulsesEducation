function pulse = PulseCreateSechSiemens(B1Max, B1Threshold, dwellTimeUs, ...
                                        VOI, gradAxis, kappa, ...
                                        eta)
% Creates a sech pulse just as I do on my Siemens ALEPH sequence.
% 
% Input Variables
% Name          Units  Description
% B1Max         kHz    Maximal B1 to which the pulse will be calibrated
% B1Threshold   kHz    The adiabatic threshold B1, above which the pulse
%                      will perform as expected.
% [dwellTimeUs] us     Default: 10 us.
% [VOI]         mm     Optional. Size of VOI (for calibrating gradient).
%                      If omitted, all gradients will be set to 0.
% [gradAxis]    'x',   Optional. The gradient axis. If omitted, the 
%               'y',   z-axis will be used.
%               'z'
% [kappa]              If omitted, set to 1.7. If present, determines the
%                      ratio between the bandwidth and the stationary
%                      region.
% [eta]                If omitted, set to 0.65. 

if (nargin<3)
    dwellTimeUs = 10;
end

if (nargin<4)
    VOI = [];
end

if (nargin<5)
    gradAxis = 'z';
end

if (nargin<6)
    kappa = 1.7;
end

if (nargin<7)
    eta = 0.65;
end

% Constants
lambda = 1.8;  % Used for calibrating the gradient
truncFactor = 5.3;

% Sweep rate at center of pulse; Also determines size of stationary region
AFPRc = (B1Threshold/eta)*(B1Threshold/eta);

% Bandwidth
AFPBW = kappa*sqrt(AFPRc);

% Duration
AFPDurationMs = truncFactor*AFPBW/AFPRc;

% Number of samples, making sure the dwell time is in multiples of
% dwellTimeUs (in us)
AFPNumSamples = floor((AFPDurationMs*1000.0/(dwellTimeUs)) + 1.0);

% Recalculate duration after rounding up the number of samples
AFPDurationMs = dwellTimeUs*0.001*AFPNumSamples;

pulse.tp = AFPDurationMs;
pulse.RFamp = zeros(1,AFPNumSamples);
pulse.RFphase = zeros(1,AFPNumSamples);
pulse.Gx = zeros(1,AFPNumSamples);
pulse.Gy = zeros(1,AFPNumSamples);
pulse.Gz = zeros(1,AFPNumSamples);

for idx=0:AFPNumSamples-1
    normTime = (idx*1.0-AFPNumSamples*0.5)/(AFPNumSamples*0.5); % From -1 to 0
    pulse.RFamp(idx+1) = 1/cosh(truncFactor*normTime);
    % pulse.RFphase(idx+1) = -mod(pi/2*(AFPDurationMs*AFPBW)/truncFactor*log(1/cosh(truncFactor*normTime)), 2*pi);
    pulse.RFphase(idx+1) = -(pi/2*(AFPDurationMs*AFPBW)/truncFactor*log(1/cosh(truncFactor*normTime)));
end

pulse.RFamp = pulse.RFamp./max(abs(pulse.RFamp(:)))*B1Max;

larmorconst = GetGyromagneticRatio('1h');
tRef  = 5.120; % ms
dxRef = 10.0; % mm
AFPRefGrad = (AFPDurationMs/tRef)*lambda*sqrt(AFPRc)/dxRef/larmorconst*1000.0; % For a 5.12 ms pulse and 10 mm slice. In mT/m
AFPDuration = round(AFPDurationMs*1000.0); % us
AFPDurationRef = 5120; 

if ~isempty(VOI)
    AFPGrad = AFPRefGrad*(AFPDurationRef/AFPDuration)*(10/VOI) * GetGyromagneticRatio('1h')/1000; % in kHz/mm
    if (nargin>3)
        switch lower(gradAxis)
            case 'x'
                pulse.Gx = ones(1,AFPNumSamples)*AFPGrad;
            case 'y'
                pulse.Gy = ones(1,AFPNumSamples)*AFPGrad;
            case 'z'
                pulse.Gz = ones(1,AFPNumSamples)*AFPGrad;
        end
    else
        pulse.Gz = ones(1,AFPNumSamples)*AFPGrad;
    end
end
