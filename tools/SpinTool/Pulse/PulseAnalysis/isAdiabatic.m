function isAdiabatic = IsAdiabatic(pulse, offset,adiabaticityFactor)
% Description: for a given offset, pulse, and REQUIRED adiabaticity factor,
% returns 1 is pulse is adiabatic throughout its duration, and 0 if not.
%
% Inputs
%
% Variables Name     Type         Units       Description
% pulse              struct       -           RF pulse structure
% offsetVec          1xN real     kHz         Offset of spin(s)
% adiabaticityFactor 1x1 real     -           Factor required of pulse
%
% Outputs
%
% isAdiabatic        1x1 boolean  -           1 if adiabatic, 0 if not [1]
%
% [1] For a pulse to be decreed adiabatic, its adiabaticity factor must
%     *exceed* the given adiabaticityFactor throughout its duration.

if (CalcMinAdiabaticityFactor(pulse,offset)<adiabaticityFactor)
    isAdiabatic = false;
else
    isAdiabatic = true;
end