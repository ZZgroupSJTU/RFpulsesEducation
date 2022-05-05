function pulse = PulseCreateAdiabatic(bandwidth, truncFactor, peakPower, duration, numSteps, pulseType)
% Description: this function generates an adiabatic 
% hyperbolic secant inversion pulse.
%
% Variables Name  Units  Description
% bandwidth       kHz    Sweep bandwidth
% truncFactor     -      Truncation factor[1]
% peakPower       kHz    Peak power of pulse
% duration        ms     Pulse total duration
% numSteps        -      Total steps
% pulseType       -      Pulse family (a string).
%                        Available choices:
%                        'sech', 'sech8'
%
% [1] The truncation factor determines when the 
%     hyperbolic secant gets cut off. Higher values 
%     mean more of the sech fits inside the pulse 
%     duration, but at the cost of "wasting" time. 
%     A reasonable value would be 5.2, since 
%     sech(5.2)=0.01


normTimeAxis = 2*linspace(0,1, numSteps) - 1;   % Ranges from -1 to 1
normDwellTime = normTimeAxis(2) - normTimeAxis(1);
dwellTime = duration/numSteps;
pulse.tp = duration;

switch (lower(pulseType))
    case 'sech'
        pulse.RFamp = peakPower*sech(truncFactor*normTimeAxis);
        pulse.RFphase = (pi/2)*(duration*bandwidth)/truncFactor*log(sech(truncFactor*normTimeAxis));
    case 'sech8'
        pulse.RFamp = peakPower*sech(truncFactor*normTimeAxis.^8);
        wRFNormalized = cumsum(sech(truncFactor*normTimeAxis.^8).^2*normDwellTime);
        wRFNormalized = wRFNormalized./(max(abs(wRFNormalized)));
        wRFNormalized = wRFNormalized - 1/2;
        wRF = 2*pi*bandwidth*wRFNormalized;
        RFphase = cumsum(wRF*dwellTime);
        RFphase = RFphase - RFphase(1);
        pulse.RFphase = RFphase;
end

pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);
