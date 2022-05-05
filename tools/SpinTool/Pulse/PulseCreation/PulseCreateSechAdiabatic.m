function pulse = PulseCreateSechAdiabatic(thresholdB1, truncFactor, peakPower, numSteps)
% Description: this function generates an adiabatic 
% hyperbolic secant inversion pulse, which inverts all spins wherever
% its amplitude meets or exceeds thresholdB1.
%
% Variables Name  Units  Description
% thresholdB1     kHz    Minimum B1 for inversion.
% truncFactor     -      Truncation factor [1, 2]
% numSteps        -      Total steps
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


BW = 2*thresholdB1; % kHz
goodnessFactor = 3;
eta = 1;
duration = goodnessFactor*2*truncFactor/(pi*eta*eta*BW); % ms

pulse = PulseCreateSech(BW, truncFactor, peakPower, duration, numSteps);