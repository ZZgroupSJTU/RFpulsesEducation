function pulse = PulseCreateAdiabaticSech(bandwidth, truncFactor, peakPower, adiabaticityFactor, spectralWidth)
% Creates an adiabatic hyperbolic secant pulse with a specified
% adiabaticity factor (which must be >> 1 to ensure adiabaticity). 
% 
% Variables Name      Units  Description
% bandwidth           kHz    Sweep bandwidth
% truncFactor         -      Truncation factor [1, 2]
% peakPower           kHz    Peak power of pulse 
% adiabaticityFactor  ms     Pulse total duration [3]
% spectralWidth       kHz    Total spectral width of sample. [4]
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
% [3] Adiabaticity is ensured by making the time*bandwidth product
%     "large enough". The meaning of "large enough" depends on the
%     ratio r = 2*peakPower/bandwidth:
%     r>1:   (pi/(2*b))*BW*T/r   = adiabaticityFactor >> 1
%     r<1:   (pi/(2*b))*BW*T*r^2 = adiabaticityFactor >> 1
%     where BW is the bandwidth in kHz, T is the duration in ms, and
%     b=5.2 is the truncation factor.
% [4] Used to set the step size (to avoid aliasing), according to
%     1/spectralWidth = dt 


r = 2*peakPower/bandwidth;
A = pi/(2*truncFactor);

if (r>=1)
    duration = r*adiabaticityFactor/(A*bandwidth);
else
    duration = adiabaticityFactor/(A*bandwidth*r^2);
end

dwellTime = 1/spectralWidth;
numSteps = ceil(duration/dwellTime);

pulse = PulseCreateSech(bandwidth, truncFactor, peakPower, duration, numSteps);
