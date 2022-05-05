function [timeAxis, adiabaticityFactor, effField, rateChange, effFieldX, effFieldZ, subtendedAngle] = CalcAdiabaticityFactor(pulse, offset)
% Function name: CalcAdiabaticityFactor
% SYNTAX: [adiabaticityFactor, effField, rateChange] = CalcAdiabaticityFactor(pulse, offset)
%
% Description: Calculates the adiabaticity factor for a given pulse &
% offset, as a function of time for the pulse's duration. An adiabaticity 
% factor much larger than 1 indicates the pulse is adiabatic. 
%
% Inputs
%
% Variables Name     Type         Units    Description
% pulse              struct       -        RF pulse structure
% offset             1x1 real     kHz      Offset of spin(s)
%
% Outputs
%
% Variables Name     Type         Units    Description
%
% rateChange         1xN real     kHz      Rate of change of FM field [1]
% effField           1xN real     kHz      Size of FM field
% adiabaticityFactor 1xN real     -        Ratio of effField/rateChange [2]
% effFieldX, 
% effFieldZ          1xN real     kHz      X and Z components of field in
%                                          FM frame
% subtendedAngle     1xN real     radians  Angle made with z-axis by the
%                                          effective field in the FM frame
%
% [1] The angle subtended by the field with the z-axis in the FM frame is
%     given by tan(psi) = (x-Component of field/z-Component of field). 
%     rateChange computes 1/(2*pi)*dpsi/dt.
% [2] For adiabaticity, this factor should be well over 1 throughout the
%     pulse.
% [3] All output quantities are given as a function of time (N = number of
%     time steps in pulse)

numSteps = length(pulse.RFamp);
dt = pulse.tp/numSteps;

timeAxis = [0:dt:(numSteps-2)*dt];

% Differentiate the RF phase to compute the instantaneous frequency in kHz
pulseFreq = 1/(2*pi)*(diff(pulse.RFphase))/dt;

effFieldZ = pulseFreq;
effFieldX = pulse.RFamp(1:end-1);

% Calculate the frequency modulated field's strength as func. of time, in kHz
effField = sqrt(pulse.RFamp(1:end-1).^2 + (offset - pulseFreq).^2);

% Calculate the rate of change of the angle the field subtends with the
% z-axis, in kHz (1/(2*pi)*dpsi/dt)
subtendedAngle = atan2(offset - pulseFreq, pulse.RFamp(1:end-1));
rateChange = 1/(2*pi)*diff(subtendedAngle)/dt;
rateChange(end+1)=rateChange(end);


adiabaticityFactor = effField./rateChange;