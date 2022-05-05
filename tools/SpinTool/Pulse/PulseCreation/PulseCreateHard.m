function pulse = PulseCreateHard(duration, tiltAngle, pulsePhase, numSteps)
% pulse = PulseCreateHard(duration, tiltAngle, pulsePhase, numSteps)
%
% The current function creates a pulse structure for a hard pulse, 
% about a given axis. The pulse's power is determined
% directly from the input parameters, according to the relation:
% gamma*[RF Power in kHz]*[Pulse Duration] = [Tilt Angle]
%
% Input:
%
% Variable Name    Units         Description
% duration         microseconds  Pulse duration
% tiltAngle        degrees       Tilt angle, in degrees
% pulsePhase       degrees       Phase of excitation axis (0 for x)
% numSteps         -             Optional. Creates the appropriate number 
%                                of steps in the pulse. If unspecified,
%                                a single step will be used.

if (nargin<4)
    numSteps = 1;
end

% Convert duration to ms
duration = duration*0.001;

% Convert tiltAngle to radians
tiltAngle = tiltAngle/360*2*pi;
pulsePhase = pulsePhase/360*2*pi;

% Create time vector
pulse.tp = duration;
B1kHz = tiltAngle/(2*pi*duration);
pulse.RFamp = ones(1,numSteps)*B1kHz;
pulse.RFphase = pulsePhase*ones(1,numSteps);

pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);