% ========================================================================
% Unit Test: Pulse-Acquire for GABA using ApplySequenceJ
% ========================================================================
%
% Simulates a "full" MEGA-PRESS experiment on a simple GABA (J-coupled)
% system, using the ApplySequenceJ routine. The MEGA-PRESS sequence is 
% (following Waddell, MRI 25:1032-1038 (2007))
%
%
% FS - frequency selective editing pulse 
%
%             180     180(FS)     180      180(FS)
%      90      _        _          _          _
%       _     | |      | |        | |        | |
%      | |    | |      | |        | |        | |
% RF __| |____| |______| |________| |________| |_______[Acquire]
%       <--t1--><--t2---><---t3----><---t4----><---t5-->
%
% Following Waddell 2007, we put
%     t1 = shortest
%     t2 = TE/4
%     t3 = TE/2-t4
%     t4 = TE/4
%     t5 = TE/2-t1-t2
% such that
%     t1+t2+t3+t4+t5     = TE
%     t1-t2+t3-t4+t5     = 0   (ON)
%     t1-(t2+t3)+(t4+t5) = 0   (OFF)
% 
% The editing and position selection pulses are full, finite time pulses.
% However, only a single spin at the center of the "voxel" (0,0,0) is
% used, so the finiteness of the pulses has negligible effect on the
% result compared to hard pulses (the comparison IS carried out). 
% The full pulses are mainly incorporated to see we can actually fulfill
% the timing constraints.
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
linewidth = 4; % Hz

spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 

spins = SpinsJAddMolecule(spins, 'gaba');
% spins = SpinsJAddMolecule(spins, 'glu');
% spins = SpinsJAddMolecule(spins, 'naa singlet');
T2 = 1000/(linewidth*pi); % Effective T2, in ms (for ad-hoc acquisition line broadening)


% ========================================================================
% Load the pulses.
% B1 maximal amplitude is calibrated to 2.5 kHz, say for a preclinical
% scanner. 
% ========================================================================

VOI = [10 10 10]; % Doesn't really matter
flipAngle = [90 180 180]; % PRESS
PRESS90Duration = 1.5; % ms
PRESS180Duration = 3; % ms

% Volume selection PRESS pulses: 90-180-180 along the x, y and z axes.
pX = PulseReadSiemensInclude('SLR90.h',  PRESS90Duration,  flipAngle(1), VOI(1), 'x');
pY = PulseReadSiemensInclude('SLR180.h', PRESS180Duration, flipAngle(2), VOI(2), 'y');
pZ = PulseReadSiemensInclude('SLR180.h', PRESS180Duration, flipAngle(3), VOI(3), 'z');

% Frequency selecting refocusing pulses (for the ON scans)
tiltAngle = 180;
SW = 4;
widthFreq = 0.06; % kHz
refPPM = 4.7; % ppm (water)
centerFreqPPM = 1.9; % ppm (GABA irradiation)
centerFreq = (centerFreqPPM - refPPM)*(B0*GetGyromagneticRatio('1h'))*0.001; % kHz
pFS = PulseCreateGaussian(centerFreq, widthFreq, SW, tiltAngle);

% ========================================================================
% Define editing sequence: Ideal & Full simulations
% ========================================================================

Gx = 0;
Gy = 0;
Gz = 0;
SW = 1.2; % kHz
numAcqPoints = 1000; 
dt = 1/SW;
timeAxis = [0:dt:(numAcqPoints-1)*dt];
TE = 68; % ms

% Calculate the inter-pulse delays, from center to center
t1 = 2.5; % Choose the "shortest possible" t1 & t2
t2 = TE/4;
t4 = TE/4;
t3 = TE/2-t4;
t5 = TE/2-t1-t2;

% Delays from end of one pulse to the beginning of the next
d1 = t1 - pX.tp/2  - pY.tp/2;
d2 = t2            - pY.tp/2  - pFS.tp/2;
d3 = t3                       - pFS.tp/2 - pZ.tp/2;
d4 = t4                                  - pZ.tp/2  - pFS.tp/2;
d5 = t5                                             - pFS.tp/2;

% Check to see delays are non-negative
if t1<0, error('Negative t1 = %.2f', t1); end
if t2<0, error('Negative t2 = %.2f', t2); end
if t3<0, error('Negative t3 = %.2f', t3); end
if t4<0, error('Negative t4 = %.2f', t4); end
if t5<0, error('Negative t5 = %.2f', t5); end

if d1<0, error('Negative d1 = %.2f', d1); end
if d2<0, error('Negative d2 = %.2f', d2); end
if d3<0, error('Negative d3 = %.2f', d3); end
if d4<0, error('Negative d4 = %.2f', d4); end
if d5<0, error('Negative d5 = %.2f', d5); end

affectedNucleiExc = {[1 1 1 1 1 1], [1 1 1 1 1], [0]};
affectedNuclei = {[1 1 0 0 0 0], [1 1 1 1 0], [1]};

% Ideal MEGA-PRESS, for comparison
seqONIdeal  = {{'hard', 90, 270, affectedNucleiExc}, {'delay', t1}, ...
               {'hard', 180, 0},                     {'delay', t2}, ...
               {'hard', 180, 0, affectedNuclei},     {'delay', t3}, ...    % Editing pulse
               {'hard', 180, 0},                     {'delay', t4}, ...
               {'hard', 180, 0, affectedNuclei},     {'delay', t5}, ...    % Editing pulse
              };
seqOFFIdeal = {{'hard', 90, 270},                    {'delay', t1}, ...
               {'hard', 180, 0},                     {'delay', t2+t3}, ...
               {'hard', 180, 0},                     {'delay', t4+t5}, ...
              };

% Full MEGA-PRESS simulation, using finite pulses
seqON       = {{'pulse', pX, [], 270},               {'delay', d1}, ...
               {'pulse', pY},                        {'delay', d2}, ...
               {'pulse', pFS},                       {'delay', d3}, ...    % Editing pulse
               {'pulse', pZ},                        {'delay', d4}, ...
               {'pulse', pFS},                       {'delay', d5}, ...    % Editing pulse
              };
seqOFF      = {{'pulse', pX, [], 270},               {'delay', d1}, ...
               {'pulse', pY},                        {'delay', d2+pFS.tp+d3}, ...
               {'pulse', pZ},                        {'delay', d4+pFS.tp+d5}, ...
              };

fprintf('Preparation report:\n');
fprintf('    Pulse durations (ms): \n');
fprintf('        PRESS: 90, 180, 180  = %.2f, %.2f, %.2f\n', pX.tp, pY.tp, pZ.tp);
fprintf('        MEGA:  180           = %.2f\n', pFS.tp);
fprintf('    Pulse amplitudes (kHz): \n');
fprintf('        PRESS: 90, 180, 180  = %.2f, %.2f, %.2f\n', max(pX.RFamp), max(pY.RFamp), max(pZ.RFamp));
fprintf('        MEGA:  180           = %.2f\n', max(pFS.RFamp));
fprintf('    Pulse gradients (mT/m): \n');
fprintf('        PRESS: 90, 180, 180  = %.2f, %.2f, %.2f\n', max(pX.Gx)*1000/42.576, max(pY.Gy)*1000/42.576, max(pZ.Gz)*1000/42.576);
fprintf('    Pulse delays (center to center, ms): \n');
fprintf('        t1, t2, t3, t4, t5   = %.2f, %.2f, %.2f, %.2f, %.2f\n', t1, t2, t3, t4, t5);
fprintf('    Pulse delays (end to start, ms): \n');
fprintf('        d1, d2, d3, d4, d5   = %.2f, %.2f, %.2f, %.2f, %.2f\n', d1, d2, d3, d4, d5);
fprintf('    Total sequence times, start of first pulse to acquisition (ms):\n');
fprintf('        Ideal, ON:  %.2f\n', CalcSeqTotalTime(seqONIdeal));
fprintf('        Ideal, OFF: %.2f\n', CalcSeqTotalTime(seqOFFIdeal));
fprintf('        Full,  ON:  %.2f\n', CalcSeqTotalTime(seqON));
fprintf('        Full,  OFF: %.2f\n', CalcSeqTotalTime(seqON));
fprintf('        Full,  ON:  %.2f (without first half of excitation pulse; i.e. this is TE)\n', CalcSeqTotalTime(seqON) - pX.tp/2);
fprintf('        Full,  OFF: %.2f (without first half of excitation pulse; i.e. this is TE)\n', CalcSeqTotalTime(seqOFF) - pX.tp/2);

% ========================================================================
% Apply sequence
% ========================================================================

tic

fprintf('Full (non-ideal) simulation ... \n');
fprintf('    Commencing ON simulation\n');
[spinsON, ~] = ApplySequenceJ(spins, seqON);

fprintf('    Commencing OFF simulation\n');
[spinsOFF, ~] = ApplySequenceJ(spins, seqOFF);

fprintf('    Creating transition tables\n');
TTON = CreateTransitionTableJ(spinsON);
TTOFF = CreateTransitionTableJ(spinsOFF);

fprintf('Ideal simulation ... \n');
fprintf('    Commencing ON simulation\n');
[spinsONIdeal, ~] = ApplySequenceJ(spins, seqONIdeal);

fprintf('    Commencing OFF simulation\n');
[spinsOFFIdeal, ~] = ApplySequenceJ(spins, seqOFFIdeal);

fprintf('    Creating transition tables\n');
TTONIdeal = CreateTransitionTableJ(spinsONIdeal);
TTOFFIdeal = CreateTransitionTableJ(spinsOFFIdeal);

fprintf('Total simulation time: %.2f sec\n', toc);

% ========================================================================
% Create FIDs from the transition tables
% ========================================================================


% Use the transition table to create the FIDs for the full MEGA-PRESS simulation
fidON = zeros(1, numAcqPoints);
fidOFF = zeros(1, numAcqPoints);

numLines = size(TTON,1);
for idx=1:numLines
    fidON = fidON + TTON(idx,2)*exp(-timeAxis/T2).*exp(-1i*2*pi*TTON(idx,1)*timeAxis);
end
numLines = size(TTOFF,1);
for idx=1:numLines
    fidOFF = fidOFF + TTOFF(idx,2)*exp(-timeAxis/T2).*exp(-1i*2*pi*TTOFF(idx,1)*timeAxis);
end
fid = {fidON, fidOFF};


% Create FID for the ideal acquisition (for comparison)
fidONIdeal = zeros(1, numAcqPoints);
fidOFFIdeal = zeros(1, numAcqPoints);

numLinesIdeal = size(TTONIdeal,1);
for idx=1:numLinesIdeal
    fidONIdeal = fidONIdeal + TTONIdeal(idx,2)*exp(-timeAxis/T2).*exp(-1i*2*pi*TTONIdeal(idx,1)*timeAxis);
end
numLinesIdeal = size(TTOFFIdeal,1);
for idx=1:numLinesIdeal
    fidOFFIdeal = fidOFFIdeal + TTOFFIdeal(idx,2)*exp(-timeAxis/T2).*exp(-1i*2*pi*TTOFFIdeal(idx,1)*timeAxis);
end
fidIdeal = {fidONIdeal, fidOFFIdeal};

fidAll = [fid, fidIdeal];

% ========================================================================
% Plot results
% ========================================================================

spec2=SPT1DSpectrum('B0', B0, 'FID', fidAll, 'dwellTime', 0.001/SW, 'title', {'ON', 'OFF', 'ON (Ideal)', 'OFF (Ideal)'});
