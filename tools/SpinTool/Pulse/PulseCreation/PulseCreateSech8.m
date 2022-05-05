function pulse = PulseCreateSech8(bandwidth, truncFactor, peakPower, duration, numSteps)
% Description: this function generates a hyperbolic secant-8 pulse.
%
% Variables Name  Units  Description
% bandwidth       kHz    Sweep bandwidth
% truncFactor     -      Truncation factor [1, 2]
% peakPower       kHz    Peak power of pulse [3]
% duration        ms     Pulse total duration
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


normTimeAxis = 2*linspace(0,1, numSteps) - 1;   % Ranges from -1 to 1
dwellTime = duration/numSteps;

pulse.tp = duration;
pulse.RFamp = peakPower*sech(truncFactor*normTimeAxis.^8);

F2 = cumsum(pulse.RFamp.^2);
F2 = F2-mean(F2);
F2 = F2./max(abs(F2));

wRF = 2*pi*(bandwidth/2)*F2;


pulse.RFphase = cumsum(wRF*dwellTime);
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);