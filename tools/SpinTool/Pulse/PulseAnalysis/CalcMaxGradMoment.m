function maxMoment = CalcMaxGradMoment(totalTime, maxGradAmp, maxSlewRate)
% Assuming a trapezoidal gradient, returns the maximal moment (in m^(-1))
% realizable in a given time interval totalTime.
%
% Input Variables
% Name           Units    Description
% totalTime      ms       Total gradient duration
% [maxGradAmp]   mT/m     Maximal gradient amplitude. Defaults to 15 mT/m.
% [maxSlewRate]  mT/m/ms  Maximal slew rate. Defaults to 100 mT/m/ms.

if nargin<2, maxGradAmp = 15; end
if nargin<3, maxSlewRate = 100; end

% % Maximum spoiling moment, given that 
% maxTriangleMoment = maxGradAmp^2/maxSlewRate;

% If we drive the gradients at max. slew rate until the maximal amplitude
% and then back down again (i.e., a triangular waveform), how long will
% it take?
maxTriangleDuration = maxGradAmp/maxSlewRate;  % ms

if (totalTime>maxTriangleDuration) % Trapezoidal gradient shape
    maxTriangleMoment = maxTriangleDuration*maxGradAmp*GetGyromagneticRatio('1h')/2; % 1/m
    maxFlatMoment = (totalTime - maxTriangleDuration)*maxGradAmp*GetGyromagneticRatio('1h'); % 1/m
    maxMoment = maxTriangleMoment + maxFlatMoment;
else % Triangular gradient shape
    maxAttainableGrad = (totalTime/2)*maxSlewRate*GetGyromagneticRatio('1h'); % kHz/m
    maxMoment = maxAttainableGrad*totalTime/2; % 1/m 
end