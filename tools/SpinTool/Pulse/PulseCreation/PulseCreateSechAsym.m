function pulse = PulseCreateSechAsym(bandwidth, truncFactor, peakPower, duration, numSteps, asymFactor)
% Description: this function generates an adiabatic 
% hyperbolic secant inversion pulse.
%
% Variables Name  Units  Description
% bandwidth       kHz    Sweep bandwidth
% truncFactor     -      Truncation factor [1, 2]
% peakPower       kHz    Peak power of pulse [3]
% duration        ms     Pulse total duration
% numSteps        -      Total steps
% asymFactor      -      Asymmetry factor [4]
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
% [4] The asymmetry factor is a number between 0 to 1, indicating 
%     at what point should the pulse shape be truncated. 1 indicates
%     no truncating. 0.5 indicates truncating halfway.


normTimeAxis = linspace(-1,-1 + asymFactor*2, numSteps);   % Ranges from -1 to 1
pulse.tp = duration;
pulse.RFamp = peakPower*sech(truncFactor*normTimeAxis);
pulse.RFphase = (pi/2)*(duration*bandwidth)/truncFactor*log(sech(truncFactor*normTimeAxis));
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);