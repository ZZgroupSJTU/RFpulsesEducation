close all
clear all
clc


% ========================================================================
% Define profile
% ========================================================================

filterOrder = 256;
freqBand = [0 0.1 0.18 1.0];
freqResponse = [1 1 0 0];
B=firls(filterOrder,freqBand,freqResponse);

%plot(abs(fftshift(fft(fftshift(B)))))

A=b2a(B);
%SLRPlotMagnetizationFromAB(A, B, [0; 0; 1], 128);

% ========================================================================
% ISLR
% ========================================================================

pulseDuration = 10;
pulse = SLRInverseTransform(A, B, pulseDuration);

pulse2 = ab2rf(A, B);

figure
subplot(2,1,1);
plot(pulse.RFamp);
subplot(2,1,2);
plot(abs(pulse2));

ISLR2 = RFToPulse(pulse2, 40);

% ========================================================================
% Simulate
% ========================================================================

numFreqPoints = 1000;
% initMag = [1; 0; 0];
initMag = [0; 0; 1];
[Mx, My, Mz, phaseMxy, magMxy, freqAxis] = CalcPulseFreqResponseWithRelax(ISLR2, 1e6, 1e6, -4, 4, numFreqPoints, initMag, 0, 0);
figure
subplot(2,2,1);
plot(freqAxis, Mz);
title('Mz');
subplot(2,2,2);
plot(freqAxis,magMxy);
title('Mxy');

