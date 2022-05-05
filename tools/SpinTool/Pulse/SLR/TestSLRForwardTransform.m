close all
clear all
clc


% ========================================================================
% Create a 90-deg. SINC pulse and calculate its frequency response profile
% ========================================================================

pulse = PulseCreateSinc(6, 5, 500);
pulse.RFamp = pulse.RFamp*0.6;
numFreqPoints = 200;
[Mx, My, Mz, phaseMxy, magMxy, freqAxis] = CalcPulseFreqResponseWithRelax(pulse, 1e6, 1e6, -8, 8, numFreqPoints, [0; 0; 1], 0, 0);


% ========================================================================
% Calculate the FSLR
% ========================================================================

[A,B] = SLRForwardTransform(pulse);

SLRPlotMagnetizationFromAB(A, B, [0; 0; 1], 512);


return


% ========================================================================
% Calculate the ISLR
% ========================================================================

pulseISLR = SLRInverseTransform(A, B, pulse.tp);

figure
% plot(abs(pulse.RFamp));
% hold
% plot(abs(pulseISLR.RFamp)./max(abs(pulseISLR.RFamp))*max(abs(pulse.RFamp)),'r');
%plot(pulse.RFamp.*cos(pulse.RFphase));
%hold
%plot(abs(pulseISLR.RFamp)./max(abs(pulseISLR.RFamp))*max(abs(pulse.RFamp)).*cos(pulseISLR.RFphase),'r');
subplot(2,1,1);
plot(real(pulse.RFamp.*exp(1i*pulse.RFphase)));
hold
pulseSLR = pulseISLR.RFamp.*exp(1i*pulseISLR.RFphase);
pulseSLR = pulseSLR(2:501);
plot(real(pulseSLR)+0.001,'r');
title('Real(B1): Original(blue), ISLR (red) ');

subplot(2,1,2);
plot(imag(pulse.RFamp.*exp(1i*pulse.RFphase)));
hold
plot(imag(pulseSLR),'r');
title('Imag(B1): Original(blue), ISLR (red) ');

