function spins = ApplyPulseRectRelax(spinsIn,pulseDuration, tiltAngle, pulsePhase)
% Applies a rectangular pulse to the spins, in the presence of relaxation.
% The rectangular pulse can be made as "hard" as needed by simply taking
% the duration to be shorter and shorter. Note that this isn't going to
% limit the peak RF power in any way.
%
% Input Parameters
% Variable Name      Units      Description
% spinsIn            -          Input spin structure
% pulseDuration      ms         Duration of rectangular pulse
% tiltAngle          degrees    Tilt angle for spins on resonance (no
%                               offset)
% pulsePhase         degrees    Phase of pulse. Take -y=270 for Mz-->Mx.
%
%
% Output Parameters
% Variable Name      Units      Description
% spinsOut           -          Output spin structure

% Create rectangular pulse
pulse = PulseCreateHard(pulseDuration*1000, tiltAngle, pulsePhase);

% Apply it
spins = ApplyPulseRelax(spinsIn, pulse);