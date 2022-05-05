function pulse = PulseCreateHSn(n, bandwidth, truncFactor, peakPower, duration, numSteps, isCorrectBase)
% Description: this function generates an adiabatic 
% hyperbolic HSn secant inversion pulse. See Garwood's paper,
% "The Return of the Adiabatic Sweep", for more details on these
% pulses (basically they are a more "squished" Sech pulse, with more SAR
% and higher bandwidth.)
%
% Variables Name  Units  Description
% n               -      Order of Sech pulse.
% bandwidth       kHz    Sweep bandwidth
% truncFactor     -      Truncation factor [1, 2]
% peakPower       kHz    Peak power of pulse [3]
% duration        ms     Pulse total duration
% numSteps        -      Total steps
% isCorrectBase   0,1    1 by default.
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

if (nargin<7), isCorrectBase = 1; end;

normTimeAxis = 2*linspace(0,1, numSteps) - 1;   % Ranges from -1 to 1
dwellTime = duration/numSteps;

pulse.tp = duration;
pulse.RFamp = peakPower*sech(truncFactor*normTimeAxis.^n);


if (isCorrectBase==1), pulse.RFamp = pulse.RFamp - pulse.RFamp(1); end;

% Calculate frequency sweep function (wRF(t))
dtau = normTimeAxis(2)-normTimeAxis(1);
F2 = cumsum(sech(truncFactor*normTimeAxis.^n).^2*dtau);
F2 = F2-mean(F2);
F2 = F2./max(abs(F2));
wRF = 2*pi*(bandwidth/2)*F2;


pulse.RFphase = -cumsum(wRF*dwellTime); % The minus makes it conform to the traditional Sech definitions in Hoult's paper
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);

