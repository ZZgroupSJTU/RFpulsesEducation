% ========================================================================
% Unit Test: A simple spin-echo using ApplySequenceJ
% ========================================================================
%
% Simulates a spin-echo experiment on a simple 2 spin (J-coupled)
% system and a simple non-coupled spin-1/2, using the ApplySequenceJ 
% routine.
%
% ========================================================================

clear all
close all
clc

% ========================================================================
% Define spin system
% ========================================================================

csCenter = 4.7; % ppm
B0 = 3; % Tesla
B1 = 1; % Scaling
isSecular = 1;
linewidth = 0; % Hz

spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 
spins = SpinsJAddMolecule(spins, 'spin half');
spins = SpinsJAddMolecule(spins, 'simple J');

% ========================================================================
% Define sequence
% ========================================================================

Gx = 0;
Gy = 0;
Gz = 0;
SW = 2; % kHz
numAcqPoints = 2000; 

TE = 100; % ms

seq = {{'hard', 90, 270}, {'delay', TE/2}, {'hard', 180, 0}, {'delay', TE/2}}; 

% ========================================================================
% Apply sequence
% ========================================================================

tic
[spins, ~] = ApplySequenceJ(spins, seq);
% Acquisition time axis
dwellTime = 1/SW; % ms
acqTimeAxis = [0:dwellTime:(numAcqPoints-1)*dwellTime]; % ms
FWHM = 2; % post hoc linewidth of each transition (in Hz)
[~, fid] = CreateTransitionTableJ(spins, 'AcqType', 'separate', 'FWHM', FWHM, 'TimeAxis', acqTimeAxis, 'Domain', 'time');
fprintf('Total time simulating: %.2f sec\n', toc)

% ========================================================================
% Process and plot results
% ========================================================================

SPT1DGUI('fid', fid, 'B0', B0, 'dwellTime', dwellTime*0.001);
