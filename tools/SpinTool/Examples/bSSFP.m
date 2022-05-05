% =========================================================================
% 
% bSSFP
%
% =========================================================================

clear all
close all
clc

% =========================================================================
% Create an ensemble of spins with a Gaussian distribution of frequencies
% =========================================================================

FWHM = 1; % Hz
numChemShifts = 1000;  
csVec = randn(1, numChemShifts)*FWHM*0.001; % kHz
T1 = 1e3; % ms
T2 = 300; % ms
M0 = 1; % a.u.
sampleSize = 10; % Doesn't matter
initMag = [0;0;1]; % Initial magnetization is at thermal eq. pointing along z
numSpinsPerShift = 1;
spins = InitSpinsRelax(csVec, numSpinsPerShift, sampleSize, initMag, T1, T2, M0);

% =========================================================================
% Create bSSFP sequence
% =========================================================================

FA = 60;  % degrees
TR = 80; % ms
numReps = 34;

dwellTime = 1; % ms
SW = 1/dwellTime;
numAcqPts = round(TR/dwellTime);
seq = repmat({{'hard', FA}, {'acquire', numAcqPts, SW, 0, 0, 0}}, [1, numReps]);


% =========================================================================
% Run the simulation
% =========================================================================

[spins, fidCellArray] = ApplySequence(spins, seq);


% =========================================================================
% Plot the resulting FID as function of time
% =========================================================================

figure
fid = cell2mat(fidCellArray);
plot(abs(fid))