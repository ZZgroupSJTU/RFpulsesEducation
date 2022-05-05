% ========================================================================
% Unit Test: Pulse-Acquire using ApplySequenceJ
% ========================================================================
%
% Simulates a pulse-acquire experiment on a simple 2 spin (J-coupled)
% system, using the ApplySequenceJ routine and possibly the 
% CreateTransitionTableJ command.
%
% ========================================================================

clear all
close all
clc

% ========================================================================
% Define spin system
% ========================================================================

% Receiver/transmitter offset, in ppm. This gets subtracted from the cs
% of each spin, in ppm, so setting this to 4.7 would, e.g., make water
% resonate at 0 Hz
csCenter = 4.7; 

% Overall B0 field, in Tesla
B0 = 3;

% Global B1 scaling. Multiplies ALL RF pulses for all spins in all molecules,
% in addition to individual B1 scaling.
B1 = 1; 

% Set to 1 to enable the secular J-coupling Hamiltonian. This slightly
% speeds up the simulation, but in general you should keep it turned off
% to simulate the spin system fully.
isSecular = 0;

% Global linewidth, in Hz. The linewidth of each molecule is then 
% calculated using:
%   Rtot = R(global) (Hz)    + R(individual) (Hz)
%        = spins.linewidth   + spins.molecule(i).spin(j).linewidth
% Note that
%   T2   = 1/(pi*Rtot)
linewidth = 1; % Hz 

% Create the spin system. This does not populate the system with 
% specific molecules or spins
spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 

% Adds molecules.
% You can change the molecules easily to anything that's allowed by 
% SpinsJAddMolecule, such as 'GABA', 'Glu', 'NAA', 'GSH', etc ... (try it!)
spins = SpinsJAddMolecule(spins, 'methanol');
spins = SpinsJAddMolecule(spins, 'ethanol');

% ========================================================================
% Define sequence
% ========================================================================

Gx = 0;
Gy = 0;
Gz = 0;
SW = 1; % kHz
numAcqPoints = 1000; 

% Propagate the spins up until acquisition. Acquisition will be handled
% using the "CreateTransitionTableJ" command.
seq = {{'hard', 90, 270}};

% ========================================================================
% Apply sequence
% ========================================================================

tic
% If we use transition tables, the command CreateTransitionTableJ will
% calculate the frequencies and (complex) amplitudes of each peak in
% the spectrum. It is then up to us to convolve these resonances with
% the suitable FID, which is done here.
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
