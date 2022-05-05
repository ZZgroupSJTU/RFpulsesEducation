function pulse = PulseCreateSech(bandwidth, truncFactor, peakPower, duration, numSteps, isCorrectBase)
% Description: this function generates an adiabatic 
% hyperbolic secant inversion pulse.
%
% Variables Name  Units  Description
% bandwidth       kHz    Sweep bandwidth
% truncFactor     -      Truncation factor [1, 2]
% peakPower       kHz    Peak power of pulse [3]
% duration        ms     Pulse total duration
% numSteps        -      Total steps
% isCorrectBase   0,1    Default: 1. If set to 1, the minimum value will
%                        be subtracted from the RF's waveform, eliminating
%                        jumps at the beginning (may be advisable to 
%                        set to 0 if truncFactor is low, < 5).
%
% [1] The truncation factor determines when the 
%     hyperbolic secant gets cut off. Higher values 
%     mean more of the sech fits inside the pulse 
%     duration, but at the cost of "wasting" time. 
%     A reasonable value would be 5.2, since 
%     sech(5.2)=0.01
% [2] The hyperbolic secant function defined here
%     coincides with the one defined in 
%     Hoult et al. Phys. Rev. A 31(4):2753-5,
%     provided one takes: 
%     beta = 2*truncFactor/(duration)
%     mu   = bandwidth/(2*beta) = bandwidth*duration/(4*truncFactor)
% [3] Hoult's paper quotes the following condition for
%     adiabaticity: if mu>=2, then as long as w1>=mu*beta
%     the profile is independent of w1.

if nargin<6
    isCorrectBase = 1;
end

normTimeAxis = linspace(-1, 1, numSteps);

pulse.tp = duration;
pulse.RFamp = peakPower*sech(truncFactor*normTimeAxis);
if (isCorrectBase==1), pulse.RFamp = pulse.RFamp - min(pulse.RFamp); end;

pulse.RFphase = (pi/2)*(duration*bandwidth)/truncFactor*log(sech(truncFactor*normTimeAxis));
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);

