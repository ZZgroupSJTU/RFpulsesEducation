% ========================================================================
% ApplySequenceExchange Unit Tests
% ========================================================================

clear all
close all
clc

% ========================================================================
% Unit Test 1: Simple two exchanging spins
% ========================================================================
%
% Two coupled spins with equal intensities precessing at
% 0 and 0.1 kHz. Relaxation will be set at T2=100 ms for both.
% The sequence is a pulse-acquire with a hard pulse 90(-y) exc
% This unit test will plot the resulting FID.
%
% Hint: play around with k and see that, when k>>cs, the peaks "merge"
% and when k<<cs the peaks are independent.
%
% ========================================================================

% Create spins
cs = 0.1; % kHz
chemShiftVec = [0 cs];  % kHz
k = 0; % kHz  (k=0 means no exchange)
exchangeMatrix = [-k k; k -k]; % kHz
T1 = [1e3 1e3]; % ms
T2 = [1e2 1e2]; % ms
eqMag = [1 1]; % a.u.
initialMag = [0; 0; eqMag(1); 0; 0; eqMag(2)];
spins = InitSpinsExchange(chemShiftVec, exchangeMatrix, initialMag, T1, T2, eqMag);

% Create sequence
SW = 1.0; % kHz
numPts = 1000; 
seq = {{'hard', 90, 270}, {'acquire', numPts, SW, 0, 0, 0}};

% Apply sequence
tic
[~, fid] = ApplySequenceExchange(spins, seq);
fprintf('Unit Test 1 Elapsed Time: %.2f \n', toc);

% Plot results

fid = fid{1};
figure
subplot(2,1,1)
dv = SW/numPts;
vv = [-SW/2:dv:SW/2-dv];
spec = fftshift(fft(fid));
maxSpec = max(abs(spec))*1.1;
plot(vv, real(spec));
title('Pulse-acquire: Spectrum (Real Part)');
ylim([-maxSpec*0.1 maxSpec]);
xlabel('kHz');
ylabel('a.u.');

% ========================================================================
% Unit Test 2
% ========================================================================
%
% We extend unit test #1 by applying an off-resonance constant irradiation
% to the spins at 0.1 kHz.
%
% ========================================================================

% Define saturation pulse parameters
satAngle = 1440; % deg.
satDuration = 1000; % ms
satPhase = 0; % deg.
satOffset = 0.1; % kHz
satSteps = 1;

spins = ReturnSpinsToThermalEquilibrium(spins);
seq = {{'rect', satAngle, satPhase, satDuration, satOffset, satSteps}, ...
       {'hard', 90, 270}, ...
       {'acquire', numPts, SW, 0, 0, 0}};

% Apply sequence
tic
[~, fid] = ApplySequenceExchange(spins, seq);
fprintf('Unit Test 2 Elapsed Time: %.2f \n', toc);
fid = fid{1};

% Plot results
subplot(2,1,2)
dv = SW/numPts;
vv = [-SW/2:dv:SW/2-dv];
spec = fftshift(fft(fid));
plot(vv, real(spec));
ylim([-maxSpec*0.1 maxSpec]);
title('Saturation-Pulse-Acquire: Spectrum (real part)');
xlabel('kHz');
