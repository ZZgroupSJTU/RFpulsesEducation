% ========================================================================
% Unit Test: PRESS Acquisition Pulse-Acquire for GABA using ApplySequenceJ
% ========================================================================
%
% Simulates a "full" PRESS experiment on a simple GABA (J-coupled)
% system, using the ApplySequenceJ routine. The PRESS sequence is 
%
%             180y     180z    
%      90x     _        _     
%       _     | |      | |    
%      | |    | |      | |    
% RF __| |____| |______| |________[Acquire]
%       <--t1--><--t2---><---t3--->
%
% For PRESS,
%     t1 = TE1/2
%     t2 = TE1/2 + TE2/2
%     t3 = TE2/2
% such that
%     t1+t2+t3      = TE1 + TE2 = TE
%     t1-t2+t3      = 0   
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

TE1 = 20; % ms
TE2 = 20; % ms
TE = TE1 + TE2; 

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
spins = SpinsJAddMolecule(spins, 'glu');
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
acqType = 'separate';
T2 = 50; % ms

% Calculate the inter-pulse delays, from center to center
t1 = TE1/2; 
t2 = TE1/2 + TE2/2;
t3 = TE2/2;

% Delays from end of one pulse to the beginning of the next
d1 = t1 - pX.tp/2  - pY.tp/2;
d2 = t2            - pY.tp/2  - pZ.tp/2;
d3 = t3                       - pZ.tp/2;

% Check to see delays are non-negative
if t1<0, error('Negative t1 = %.2f', t1); end
if t2<0, error('Negative t2 = %.2f', t2); end
if t3<0, error('Negative t3 = %.2f', t3); end

if d1<0, error('Negative d1 = %.2f', d1); end
if d2<0, error('Negative d2 = %.2f', d2); end
if d3<0, error('Negative d3 = %.2f', d3); end

% Ideal PRESS, for comparison
seqIdeal = {{'hard', 90, 270},                    {'delay', t1}, ...
            {'hard', 180, 0},                     {'delay', t2}, ...
            {'hard', 180, 0},                     {'delay', t3}, ...
           };

% Full PRESS simulation, using finite pulses
seq      = {{'pulse', pX, [], 270},               {'delay', d1}, ...
            {'pulse', pY},                        {'delay', d2}, ...
            {'pulse', pZ},                        {'delay', d3}, ...
           };

fprintf('Preparation report:\n');
fprintf('    Pulse durations (ms): \n');
fprintf('        PRESS: 90, 180, 180  = %.2f, %.2f, %.2f\n', pX.tp, pY.tp, pZ.tp);
fprintf('    Pulse amplitudes (kHz): \n');
fprintf('        PRESS: 90, 180, 180  = %.2f, %.2f, %.2f\n', max(pX.RFamp), max(pY.RFamp), max(pZ.RFamp));
fprintf('    Pulse gradients (mT/m): \n');
fprintf('        PRESS: 90, 180, 180  = %.2f, %.2f, %.2f\n', max(pX.Gx)*1000/42.576, max(pY.Gy)*1000/42.576, max(pZ.Gz)*1000/42.576);
fprintf('    Pulse delays (center to center, ms): \n');
fprintf('        t1, t2, t3           = %.2f, %.2f, %.2f\n', t1, t2, t3);
fprintf('    Pulse delays (end to start, ms): \n');
fprintf('        d1, d2, d3           = %.2f, %.2f, %.2f\n', d1, d2, d3);
fprintf('    Total sequence times, start of first pulse to acquisition (ms):\n');
fprintf('        Ideal:  %.2f\n', CalcSeqTotalTime(seqIdeal));
fprintf('        Full:   %.2f\n', CalcSeqTotalTime(seq));

% ========================================================================
% Apply sequence
% ========================================================================

tic

fprintf('Full (non-ideal) simulation ... \n');
[spinsOut, ~] = ApplySequenceJ(spins, seq);

fprintf('    Creating transition tables\n');
[TT, fid] = CreateTransitionTableJ(spinsOut, acqType, timeAxis, T2);
fprintf('    Done (%.2f sec)\n', toc);

tic
fprintf('Ideal simulation ... \n');
[spinsOutIdeal, ~] = ApplySequenceJ(spins, seqIdeal);

fprintf('    Creating transition tables\n');
% [TT, fid] = CreateTransitionTableJ(spins, acqType, timeAxis, T2)
[TTIdeal, fidIdeal] = CreateTransitionTableJ(spinsOutIdeal, acqType, timeAxis, T2);
fprintf('    Done (%.2f sec)\n', toc);


fidAll = [fid, fidIdeal];

% ========================================================================
% Plot results
% ========================================================================

SPT1DGUI('B0', B0, 'FID', fidAll, 'dwellTime', 0.001/SW);
