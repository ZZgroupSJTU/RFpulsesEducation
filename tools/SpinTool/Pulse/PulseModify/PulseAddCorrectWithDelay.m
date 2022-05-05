function PulseOut = PulseAddCorrectWithDelay(pulse, centerFreq1, centerFreq2, phase1, phase2, delay)
% PulseOut = PulseAdd(pulse1,pulse2);
%
% Adds two RF pulses together. Pulses are assumed to have the same duration
% and number of steps. Gradients are set to 0. Correction is applied to
% pulse frequencies according to Steffen et. al., JMR 146:369-74 (2000).

% Before adding, the frequency-shift corrections must be added to each
% pulses (to first order), according to the JMR paper

pulse1 = PulseShiftOffset(pulse, centerFreq1);
pulse2 = PulseShiftOffset(pulse, centerFreq2);

pulse1.RFphase = pulse1.RFphase + phase1;
pulse2.RFphase = pulse2.RFphase + phase2;

dt = pulse.tp/length(pulse.RFamp);

freqCorrect12 = 2*pi*(pulse1.RFamp).^2./(2*centerFreq2);
freqCorrect21 = 2*pi*(pulse2.RFamp).^2./(2*centerFreq1);

phaseCorrect12 = -1/2*cumsum(freqCorrect12)*dt;
phaseCorrect21 = -1/2*cumsum(freqCorrect21)*dt;

pulse1.RFphase = pulse1.RFphase;% + phaseCorrect21;
pulse2.RFphase = pulse2.RFphase;% + phaseCorrect12;

% Now, add!
PulseOut = PulseAddWithDelay(pulse1, pulse2, delay);
