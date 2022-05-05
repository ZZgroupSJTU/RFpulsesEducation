% ========================================================================
% Unit Test: Does strong J-coupling affect longitudinal spin states (IzSz)?
% ========================================================================
%
% Starts from IzSz for a 2 spin system with strong J-coupling, and checks
% whether the system's density matrix evolves in any appreciable way.
%
% ========================================================================

clear all
close all
clc

% ========================================================================
% Define spin system
% ========================================================================

csCenter = 0; % ppm
B0 = 3; % Tesla
B1 = 1; % Scaling
isSecular = 0;
linewidth = 1; % Hz
J = 50; % Hz

spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 
spins = SpinsJAddMolecule(spins, 'Simple J');

spins.molecule(1).csVec = [-0.1 0.1]; 
spins.molecule(1).spin.rho = kron(Iz, I);
spins.molecule(1).JMatrix = [0 J; J 0];

fprintf('Separation between offsets (in Hz): %.2f \n', (spins.molecule(1).csVec(2) - spins.molecule(1).csVec(1))*B0*GetGyromagneticRatio('1h'));
fprintf('J Coupling: %.1f Hz\n', J);

fprintf('Density matrix of system before delay: \n');
PrintDensityMatrix(spins.molecule(1).spin.rho)

% ========================================================================
% Define sequence
% ========================================================================

delay = 40; % ms
seq = {{'delay', delay}};

% ========================================================================
% Apply sequence
% ========================================================================

[spins, fid] = ApplySequenceJ(spins, seq);

% ========================================================================
% Process and plot results
% ========================================================================

fprintf('Density matrix of system after %.0f ms delay: \n', delay);

PrintDensityMatrix(spins.molecule(1).spin.rho)