clear all
close all
clc

% Create an excitation pulse
pulse = PulseCreateSLR90(30, 1, 10, 'linearR6', 0);

% Get its frequency response
T1 = 10000;
T2 = 10000;
freqMin = -10;
freqMax = 10;
numResPoints = 1000;
initMag = [0; 0; 1];
[Mx, My, Mz, phaseMxy, magMxy, freqAxis] = CalcPulseFreqResponseWithRelax(pulse, T1, T2, freqMin, freqMax, numResPoints, initMag, 0, 0);

numSteps = length(pulse.RFamp);
dwellTime = pulse.tp/numSteps;
pulseSW = 1/dwellTime % in kHz


%plot(freqAxis, Mz);

[A,B] = SLRForwardTransform(pulse);

p = linspace(0, 2*pi, numResPoints);
z = exp(1i*p);
f = zeros(1,numResPoints);
for n=1:numResPoints
    f(n) = 0;
    for k=1:length(B)
         f(n) = f(n) + B(k)*z(n)^(-(k-1));
        %f(n) = f(n) + A(k)*z(n)^(-(k-1));
    end
end


% pulse2 = SLRInverseTransform(A,B, pulse.tp);
% pulse2.RFamp = 2*pulse2.RFamp;

% PulsePlot(pulse);
% PulsePlot(pulse2);

% [Mx2, My2, Mz2, phaseMxy2, magMxy2, freqAxis2] = CalcPulseFreqResponseWithRelax(pulse2, T1, T2, freqMin, freqMax, numResPoints, initMag, 0, 0);
% figure
% plot(freqAxis2, Mz2);


%figure
%subplot(3,1,1)
%plot(freqAxis,real(f));
%subplot(3,1,2)
%plot(freqAxis,imag(f));
%subplot(3,1,3)
%plot(freqAxis,abs(f));

%SLRPlotUnitCircleProfile(pulse, numResPoints);