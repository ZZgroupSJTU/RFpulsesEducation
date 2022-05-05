clear all
close all
clc

% Create a sinc pulse
N = 256;
% pulse = PulseCreateSinc(4, 12.8, N, 180);
pulse = PulseCreateSechAdiabatic(0.5, 5.2, 1.0, N);
G0 = 0.1;
pulse.Gz = G0*ones(1,N);
% Create the modulation function (don't worry about total area)
modFun = (1 - 0.5*sin([0:N-1]/(N-1)*pi));

chemShift = 0.3;

% VERSE-ify it.
pulseVERSE = PulseVERSE(pulse, modFun);
% figure
% plot(pulse.RFamp);
% hold
% plot(pulseVERSE.RFamp, 'r');
% PulsePlot({pulse, pulseVERSE})
% PulsePlot(pulseVERSE)
PlotPulseFreqResponse({pulse, pulseVERSE}, [0;0;1], -40, 40, 501, 'mz', 'z', chemShift);
