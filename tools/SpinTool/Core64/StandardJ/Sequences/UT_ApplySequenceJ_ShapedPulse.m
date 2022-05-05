% ========================================================================
% Unit Test: Shaped Pulses
% ========================================================================
%
% Simulates a shaped, frequency selective, narrow bandwidth pulse on a 
% system with many equally spaced singlets. 
%
%    180 (FS)   90
%       _        
%      | |      _ 
%      | |     | |
% RF __| |_____| |_[Acquire]
%         <-d->
%
% ========================================================================

clear all
close all
clc

% ========================================================================
% Define spin system
% ========================================================================

csCenter = 4.7; % ppm
B0 = 2.9; % Tesla
B1 = 1; % Scaling
isSecular = 0;
linewidth = 0; % Hz

spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 

csVec = 1.9;

% Add the molecule's data
LW = 4; % Linewidth, in Hz
for idx=1:numel(csVec)
    spins.molecule(idx).csVec = csVec(idx);
    spins.molecule(idx).JMatrix = 0;
    spins.molecule(idx).gmRatio = GetGyromagneticRatio('1h');
    spins.molecule(idx).spin(1).r = [0;0;0];
    spins.molecule(idx).spin(1).rho = IzN(1,1);
    spins.molecule(idx).spin(1).B1 = 1;
    spins.molecule(idx).spin(1).B0 = 0;
    spins.molecule(idx).spin(1).linewidth = LW; % Not used
end

T2 = (1/(LW*pi))*1000; % Effective T2, in ms (for ad-hoc acquisition line broadening)

% ========================================================================
% Create a narrow, frequency selective pulse
% ========================================================================

% Frequency selecting refocusing pulses (for the ON scans)
tiltAngle = 180;
pulseSW = 4;
widthFreq = 0.1; % kHz
refPPM = 4.7; % ppm (water)
centerFreqPPM = 1.9; % ppm (GABA irradiation)
centerFreq = (centerFreqPPM - refPPM)*(B0*GetGyromagneticRatio('1h'))*0.001; % kHz
pFS = PulseCreateGaussian(centerFreq, widthFreq, pulseSW, tiltAngle);

% ========================================================================
% Define editing sequence: Ideal & Full simulations
% ========================================================================

Gx = 0;
Gy = 0;
Gz = 0;
SW = 2; % kHz
numAcqPoints = 1000; 
dt = 1/SW;
timeAxis = [0:dt:(numAcqPoints-1)*dt];

% Delay between 180 FS pulse and acquisition pulse
d = 0.1; % ms

% Pulse sequence. REMEMBER THE QM SIMULATION USES THE RIGHT HAND RULE!
seq = {{'pulse', pFS}, {'delay', d}, {'hard', 90, 90}}; 

% ========================================================================
% Apply sequence
% ========================================================================

tic
[spinsOut, ~] = ApplySequenceJ(spins, seq);
TT = CreateTransitionTableJ(spinsOut);
fprintf('Total simulation time: %.2f sec\n', toc);

% ========================================================================
% Create FIDs from the transition tables
% ========================================================================

% Use the transition table to create the FIDs for the full MEGA-PRESS simulation
fid = zeros(1, numAcqPoints);
numLines = size(TT,1);
for idx=1:numLines
    fid = fid + TT(idx,2)*exp(-timeAxis/T2).*exp(-1i*2*pi*TT(idx,1)*timeAxis);
end

% ========================================================================
% Plot results
% ========================================================================

spec2=SPT1DSpectrum('B0', B0, 'FID', fid, 'dwellTime', 0.001/SW);
