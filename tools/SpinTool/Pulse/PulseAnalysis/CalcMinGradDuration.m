function [minDuration, rampTime, flatTime, gradAmp] = CalcMinGradDuration(gradMoment, maxGradAmp, maxSlewRate, nucleus)
% SYNTAX: 
%
%   minDuration = CalcMinGradDuration(gradMoment, maxGradAmp, maxSlewRate)
%
% Calculates the minimum time needed for a trapezoidal gradient, with the
% maximal gradient amplitude & slew rate given by maxGradAmp (mT/m) and 
% maxSlewRate (mT/m/ms)
%
% Input Variables
% Variable Name   Units     Description
% gradMoment      m^(-1)    Desired gradient moment. An NxM matrix, having
%                           M spoiling "periods", with N gradients along
%                           the N axes (usually 1, 2 or 3) in each period.
% [maxGradAmp]    mT/m      Max. allowed gradient amplitude. Default: 16.0
% [maxSlewRate]   mT/m/ms   Max. allowed gradient slew rate. Default: 40.0
% [nucleus]       -         '1H', '31P': determines the gyromagnetic ratio
%                           Default: '1H'
% 
% Output Variables
% Variable Name   Units     Description
% minDuration     ms        The minimum time needed to realize gradient
%                           trapezoidal shape.
% rampTime        ms        Duration of the ramp up/down time
% flatTime        ms        Duration of the flat top.
% gradAmp         mT/m      Maximal gradient reached during shape.

if (nargin<4), nucleus = '1h'; end
if (nargin<3), maxSlewRate = 40; end
if (nargin<2), maxGradAmp = 16.0; end

% The number of axes is the number of simultaneous gradient moments.
% 
numAxes = size(gradMoment, 1);
numPeriods = size(gradMoment, 2);

minDuration = zeros(1, numPeriods);
rampTime = zeros(size(gradMoment));
flatTime = zeros(size(gradMoment));
gradAmp = zeros(size(gradMoment));

gmr = GetGyromagneticRatio(nucleus);
maxRampTime = 2*maxGradAmp/maxSlewRate; % Assuming the grad. reaches maxGradAmp
maxTriangleMoment = maxGradAmp^2/maxSlewRate*gmr; % in m^(-1)

for idxPeriod=1:numPeriods
    for idxAxis=1:numAxes
        % Assumption 1: We drive the gradient's ramp up/down at full slew rate 
        % Assumtpion 2: If we have a flat top, it will be at maxGradAmp.

        % There are two possible cases: the gradient reaches the flat top
        % (trapezoidal shape) or doesn't (triangular)
        isTriangle = (maxTriangleMoment>=gradMoment(idxAxis, idxPeriod)); % 0 or 1
        if (isTriangle)
            curRampTime = sqrt(gradMoment(idxAxis, idxPeriod)/maxSlewRate/gmr);
            curFlatTime = 0;
            curGradAmp = maxSlewRate*curRampTime;
            curMinDuration = curRampTime*2;
        else
            curRampTime = maxRampTime;
            curGradAmp = maxGradAmp;
            flatMoment = gradMoment(idxAxis, idxPeriod) - maxTriangleMoment; % In m^(-1)
            curFlatTime = flatMoment/maxGradAmp/gmr;
            curMinDuration = 2*curRampTime + curFlatTime;
        end

        if curMinDuration>minDuration(idxPeriod)
            minDuration(idxPeriod) = curMinDuration;
            rampTime(idxAxis, idxPeriod) = curRampTime;
            flatTime(idxAxis, idxPeriod) = curFlatTime;
            gradAmp(idxAxis, idxPeriod) = curGradAmp;
        end
    end
end