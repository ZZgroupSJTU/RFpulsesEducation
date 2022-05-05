function pulse = PulseCreateZero(pulseDuration,numSteps)
% Description: Creates a "zero" pulse, having zero RF and gradients, with
% a given number of steps and a given duration. Such pulses can be used
% as either (i.) delays, or (ii.) in conjuction with an acquisition
% routine, which expects a pulse structure specifying the acquisition
% time points.
%
% Inputs:
%
% Variable Name   Units   Description
% pulseDuration   ms      Duration of pulse. 
% numSteps        -       Number of steps in pulse.
% 
% Outputs:
%
% Variable Name   Units   Description
% pulse           -       Output pulse structure.

pulse.tp      = pulseDuration;
pulse.RFamp   = zeros(1,numSteps);
pulse.RFphase = zeros(1,numSteps);
pulse.Gx      = zeros(1,numSteps);
pulse.Gy      = zeros(1,numSteps);
pulse.Gz      = zeros(1,numSteps);