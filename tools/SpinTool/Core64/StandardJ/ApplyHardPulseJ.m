function spins = ApplyHardPulseJ(spins, pulseDuration, tiltAngle, pulsePhase, affectedNuclei, freqRange)
% The current function applies a hard pulse to a J-coupled spin structure.
% The field's strength is calibrated according to Protons.
%
% Input:
%
% Variable Name    Units         Description
% spins            -             J-coupled spin array
% pulseDuration    ms            Pulse duration
% tiltAngle        degrees       Tilt angle, in degrees
% pulsePhase       degrees       Phase of excitation axis (0 for x)
% affectedNuclei   0, 1          Cell array of boolean vectors indicating 
%                                which spins will be affected by the pulse

if (nargin<4)
    pulsePhase = 0;
end
if nargin<5, affectedNuclei = []; end
if nargin<6, freqRange = []; end

tiltAngle = tiltAngle/180*pi; % rad
pulse.tp = pulseDuration; % ms
pulse.RFamp = tiltAngle/(2*pi)/pulse.tp; % in kHz
pulse.RFphase = pulsePhase/180*pi;
pulse.Gx = 0;
pulse.Gy = 0;
pulse.Gz = 0;
isAcquire = 0;
% PropagateJ(spins, pulse, acqType, affectedNuclei, freqRange, gradUnits, RFUnits)
[~, spins] = PropagateJ(spins, pulse, isAcquire, affectedNuclei, freqRange);
