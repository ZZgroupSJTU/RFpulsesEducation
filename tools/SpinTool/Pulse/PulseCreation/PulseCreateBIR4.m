function pulseBIR4 = PulseCreateBIR4(peakPower, tiltAngle, adiabaticityFactor, numSteps, truncFactor)
% SYNTAX:
%
%    pulseBIR4 = PulseCreateBIR4(peakPower, tiltAngle, adiabaticityFactor, numSteps, truncFactor)
%
% Creates a BIR4 adiabatic pulse of arbitrary tilt angle.

tiltAngle = tiltAngle/180*pi;

phaseBetweenAHP1 = pi + tiltAngle/2;
phaseBetweenAHP2 = - phaseBetweenAHP1;

% ========================================================================
% BIR-4 Simulation
% ========================================================================
%
% The BIR-4 can be thought of as two BIR-1 pulses, which themselves can
% be thought of as two AHP pulses:
%
%   Pulse #   4      3             2     1      
%   Phase   theta  theta/2      theta/2  0
%
%   BIR-1 =  AHP * AHP^(-1)   *   AHP * AHP^(-1)  
%            /|\   /|\            /|\   /|\
%             |_____|              |_____|
%                                                
%            Phase shift          Phase shift
%              theta/2              theta/2
%
% The overall effect of this is to produce a theta rotation about the 
% x-axis (assuming the phase of the first AHP^(-1) is 0) by an angle 
% theta.
%
% We use a sech pulse, truncated at its peak, with the following
% parameters:
%
% Amplitude (kHz)
% /|\
%  |          ___     /|\    
%  |        _/         |
%  |      _/         peakPower
%  |    _/             |
%  |  _/               |
%  |_/________________\|/_____________> time
%
%  <------------->
%   pulseDuration
%
% ========================================================================


% Compute derived parameters
bandwidth = 2*peakPower; % kHz. This is actually not the "real" bandwidth of the pulse, but of the hyperbolic secant sub-pulse 
duration = adiabaticityFactor*2*truncFactor/(pi*bandwidth); % ms
numBaseSteps = ceil(numSteps/4); % Make sure it's an integer

% Create base pulse
isCorrectBase = 0;
% pulseSech = PulseCreateSech(bandwidth, truncFactor, peakPower, duration, numBaseSteps*2, isCorrectBase);
pulseSech = PulseCreateHSn(1, bandwidth, truncFactor, peakPower, duration, numBaseSteps*2, isCorrectBase);

% Pulse 1: inverse AHP, sweeping from +x to +z in the FM frame
pulseAHP1 = pulseSech;
pulseAHP1.RFamp = pulseAHP1.RFamp(numBaseSteps+1:end);
pulseAHP1.RFphase = pulseAHP1.RFphase(numBaseSteps+1:end);
pulseAHP1.Gx = zeros(1,numBaseSteps);
pulseAHP1.Gy = zeros(1,numBaseSteps);
pulseAHP1.Gz = zeros(1,numBaseSteps);

% Pulse 2: AHP, sweeping from -z to an axis in the xy plane (in the FM 
% frame) making an angle theta/2 with the x-axis.
pulseAHP2 = pulseSech;
pulseAHP2.RFamp = pulseAHP2.RFamp(1:numBaseSteps);
pulseAHP2.RFphase = pulseAHP2.RFphase(1:numBaseSteps) + phaseBetweenAHP1;
pulseAHP2.Gx = zeros(1,numBaseSteps);
pulseAHP2.Gy = zeros(1,numBaseSteps);
pulseAHP2.Gz = zeros(1,numBaseSteps);

% Pulse 3: Inverse AHP, sweeping from axis w/ angle theta/2 in x-plane to
% the z-axis
pulseAHP3 = pulseSech;
pulseAHP3.RFamp = pulseAHP3.RFamp(numBaseSteps+1:end);
pulseAHP3.RFphase = pulseAHP3.RFphase(numBaseSteps+1:end) + phaseBetweenAHP1;
pulseAHP3.Gx = zeros(1,numBaseSteps);
pulseAHP3.Gy = zeros(1,numBaseSteps);
pulseAHP3.Gz = zeros(1,numBaseSteps);

% Pulse 4:  AHP, sweeping from +z to axis w/ angle theta in x-plane 
pulseAHP4 = pulseSech;
pulseAHP4.RFamp = pulseAHP4.RFamp(1:numBaseSteps);
pulseAHP4.RFphase = pulseAHP4.RFphase(1:numBaseSteps) + phaseBetweenAHP1 + phaseBetweenAHP2;
pulseAHP4.Gx = zeros(1,numBaseSteps);
pulseAHP4.Gy = zeros(1,numBaseSteps);
pulseAHP4.Gz = zeros(1,numBaseSteps);

% % Pulse 1: inverse AHP, sweeping from +x to +z in the FM frame
% normTimeAxis1 = linspace(0,1, numSteps);
% pulseAHP1.tp = pulseDuration;
% pulseAHP1.RFamp = peakPower*sech(truncFactor*normTimeAxis1);
% pulseAHP1.RFphase = (pi/2)*(pulseDuration*bandwidth)/truncFactor*log(sech(truncFactor*normTimeAxis1));
% pulseAHP1.Gx = zeros(1,numSteps);
% pulseAHP1.Gy = zeros(1,numSteps);
% pulseAHP1.Gz = zeros(1,numSteps);
% 
% % Pulse 2: AHP, sweeping from -z to an axis in the xy plane (in the FM 
% % frame) making an angle theta/2 with the x-axis.
% normTimeAxis2 = linspace(-1,0,numSteps);
% pulseAHP2.tp = pulseDuration;
% pulseAHP2.RFamp = peakPower*sech(truncFactor*normTimeAxis2);
% pulseAHP2.RFphase = (pi/2)*(pulseDuration*bandwidth)/truncFactor*log(sech(truncFactor*normTimeAxis2)) + phaseBetweenAHP1;
% pulseAHP2.Gx = zeros(1,numSteps);
% pulseAHP2.Gy = zeros(1,numSteps);
% pulseAHP2.Gz = zeros(1,numSteps);
% 
% % Pulse 3: Inverse AHP, sweeping from axis w/ angle theta/2 in x-plane to
% % the z-axis
% normTimeAxis3 = linspace(0,1, numSteps);
% pulseAHP3.tp = pulseDuration;
% pulseAHP3.RFamp = peakPower*sech(truncFactor*normTimeAxis3);
% pulseAHP3.RFphase = (pi/2)*(pulseDuration*bandwidth)/truncFactor*log(sech(truncFactor*normTimeAxis3)) + phaseBetweenAHP1;
% pulseAHP3.Gx = zeros(1,numSteps);
% pulseAHP3.Gy = zeros(1,numSteps);
% pulseAHP3.Gz = zeros(1,numSteps);
% 
% % Pulse 4:  AHP, sweeping from +z to axis w/ angle theta in x-plane 
% normTimeAxis4 = linspace(-1,0,numSteps);
% pulseAHP4.tp = pulseDuration;
% pulseAHP4.RFamp = peakPower*sech(truncFactor*normTimeAxis4);
% pulseAHP4.RFphase = (pi/2)*(pulseDuration*bandwidth)/truncFactor*log(sech(truncFactor*normTimeAxis4)) + phaseBetweenAHP1 + phaseBetweenAHP2;
% pulseAHP4.Gx = zeros(1,numSteps);
% pulseAHP4.Gy = zeros(1,numSteps);
% pulseAHP4.Gz = zeros(1,numSteps);


% Concatenate the 4 AHPs to form the BIR4 
pulseBIR4 = PulseConcat(pulseAHP1, pulseAHP2);
pulseBIR4 = PulseConcat(pulseBIR4, pulseAHP3);
pulseBIR4 = PulseConcat(pulseBIR4, pulseAHP4);
