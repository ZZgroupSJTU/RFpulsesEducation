function minAdiabaticityFactor = CalcMinAdiabaticityFactor(pulse, offset)
% Function name: CalcMinAdiabaticityFactor
% SYNTAX: minAdiabaticityFactor = CalcMinAdiabaticityFactor(pulse, offset)
%
% Description: Calculates the MINIMAL adiabaticity factor for a given pulse &
% offset throughout the pulse's duration.
%
% Variables Name     Type         Units       Description
% <=INPUTS=>
% pulse              struct       -           RF pulse structure
% offsetVec          1xN real     kHz         Offset of spin(s)
% <=OUPUTS=>
% minadiabaticity...
%    ...Factor       1x1 real     -           Min. (over time) of pulse's
%                                             adiabaticity factor
%
% [1] The angle subtended by the field with the z-axis in the FM frame is
%     given by tan(psi) = (x-Component of field/z-Component of field). 
%     rateChange computes 1/(2*pi)*dpsi/dt.
% [2] For adiabaticity, this factor should be well over 1.

[timeAdiabaticityFactor, effField, rateChange] = CalcAdiabaticityFactor(pulse, offset);
minAdiabaticityFactor = min(timeAdiabaticityFactor);