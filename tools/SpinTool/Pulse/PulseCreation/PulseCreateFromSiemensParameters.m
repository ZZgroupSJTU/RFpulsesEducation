function pulse = PulseCreateFromSiemensParameters(ampVec, phaseVec, ampInt, refGrad, gradAxis, pulseDuration, sliceThickness, flipAngle)
% Creates a pulse structure from an amplitude vector & phase vector, and 
% the supplied parameters, following the same calibration rules as a
% Siemens scanner (VB13 - VB17, to the best of my knowledge).
%
% Inputs:
%
% Variable Name   Units   Description
% ampVec          -       Vector containing normalized amplitudes.
% phaseVec        rad     Vector of pulse phases
% refGrad         mT/m    Reference gradient
% gradAxis        -       'slice', 'phase' or 'read'
% pulseDuration   ms 
% sliceThickness  mm 
% flipAngle       deg


numSamples = length(ampVec);

% Normalize amplitude vector, in case it's not!
ampVec = ampVec./max(abs(ampVec));

% The reference amp for the Siemens is 0.5 kHz: a 180 deg. flip in 1 ms
refAmp = 0.5;
refDuration = 5.12; % ms
refThickness = 10; % mm

% Calibrate maximal amplitude, in Khz
maxAmp = (flipAngle/180)*(numSamples/(ampInt*pulseDuration))*refAmp;

% Calibrate gradient, in mT/m
grad = refDuration/pulseDuration * (refThickness/sliceThickness) * refGrad;

% Convert gradient to kHz/mm
grad    = grad*GetGyromagneticRatio('1h')/1000;

% Create pulse structure
pulse.tp = pulseDuration;
pulse.RFamp = maxAmp.*ampVec;
pulse.RFphase = phaseVec;
pulse.Gx = zeros(1, numSamples);
pulse.Gy = zeros(1, numSamples);
pulse.Gz = zeros(1, numSamples);
switch lower(gradAxis)
    case 'z'
        pulse.Gz = ones(1, numSamples)*grad; % kHz/mm
    case 'y'
        pulse.Gy = ones(1, numSamples)*grad; % kHz/mm
    case 'x'
        pulse.Gx = ones(1, numSamples)*grad; % kHz/mm
end