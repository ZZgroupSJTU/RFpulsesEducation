function [pulse, sw, Ga, Ta] = PulseCreate_EPSI(FOV, Nz, sw_EPSI, Nv);
% function [pulse, sw, Ga, Ta] = PulseCreate_EPSI(FOV, Nz, sw_EPSI, Nv);
%
% This function generates the gradients (in a pulse) used in an EPSI experiment.
% FOV is in mm. sw_EPSI is in kHz.
%
% EPSI Illustration:
%   
%   Nv
%   /                           \
%   |     Dk           Dk       |
%   |-- Acquire ---- Acquire -- |   
%   |                           |
%   | +---------+               |
%   | | +Ga, Dt |               |
%   |-+         +--+         +- |
%   |              | -Ga, Dt |  |
%   |              +---------+  |
%   \                           /
%
% Derivation outline:
% 1. From sw  = 1/(2*Dt), compute Dt.
% 2. From dz  = FOV/Nz = 1/(Ga*Dt) compute Ga.
% 3. From FOV = 1/dk = sw/Ga compute sw ( = 1/dwell time).
% 4. Compute Ta = Nv*2*Dt

Dt = 1/(2*sw_EPSI);
Ga = Nz/(FOV*Dt);
sw = Ga*FOV;
Ta = Nv*2*Dt;

% Initialize acquisition pulse to 0
pulse.tp      = Ta;
pulse.RFamp   = zeros(1,Nv*Nz*2);
pulse.RFphase = zeros(1,Nv*Nz*2);
pulse.Gx      = zeros(1,Nv*Nz*2);
pulse.Gy      = zeros(1,Nv*Nz*2);
pulse.Gz      = zeros(1,Nv*Nz*2);

% Create & assign gradient pattern
grad = [];
for k=1:Nv
    grad = [grad,Ga*ones(1,Nz),-Ga*ones(1,Nz)];
end;

pulse.Gz = grad;
