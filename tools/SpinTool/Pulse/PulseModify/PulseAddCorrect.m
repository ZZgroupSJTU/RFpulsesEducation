function PulseOut = PulseAddCorrect(pulse, centerFreq1, centerFreq2, phase1, phase2, pulseAxis)
% SYNTAX:
%
%    PulseOut = PulseAddCorrect(pulse, centerFreq1, centerFreq2, phase1, phase2)
%
% Adds two RF pulses together. Pulses are assumed to have the same duration
% and number of steps. Gradients are set to 0. Correction is applied to
% pulse frequencies according to Steffen et. al., JMR 146:369-74 (2000).
%
% Example:
%
% % Create an SLR excitation 90-deg. pulse with peak B1 of 0.46 kHz, 
% pulse = PulseCreateSLR90(150, 0.46, 10, 'linearR18', 0);
%
% % Create two shifted copies at +1 kHz and -1 kHz
% centerFreq1 = 1; % kHz
% centerFreq2 = -1; % kHz
% pulse1 = PulseShiftOffset(pulse, centerFreq1);
% pulse2 = PulseShiftOffset(pulse, centerFreq2);
%
% % Add the two pulses without applying the correction
% pulseAddWithoutCorrection = PulseAdd(pulse1, pulse2);
% 
% % Add the two pulses with applying the correction
% pulseAddWithCorrection = PulseAddCorrect(pulse, centerFreq1, centerFreq2, 0, 0);
%
% % Plot the two on top of each other for comparison
% PlotPulseFreqResponse({pulseAddWithCorrection, pulseAddWithoutCorrection}, [0; 0; 1], -5, 5, 500, 'mz');

if nargin<6, pulseAxis = []; end

% Before adding, the frequency-shift corrections must be added to each
% pulses (to first order), according to the JMR paper

if ~isempty(pulseAxis)
    switch lower(pulseAxis)
        case 'x'
            centerFreq1 = centerFreq1*pulse.Gx(1);
            centerFreq2 = centerFreq2*pulse.Gx(1);
        case 'y'
            centerFreq1 = centerFreq1*pulse.Gy(1);
            centerFreq2 = centerFreq2*pulse.Gy(1);
        case 'z'
            centerFreq1 = centerFreq1*pulse.Gz(1);
            centerFreq2 = centerFreq2*pulse.Gz(1);
        otherwise
            % Do nothing
    end
end

pulse1 = PulseShiftOffset(pulse, centerFreq1);
pulse1.RFphase = pulse1.RFphase + phase1/180*pi;

pulse2 = PulseShiftOffset(pulse, centerFreq2);
pulse2.RFphase = pulse2.RFphase + phase2/180*pi;

dt = pulse.tp/length(pulse.RFamp);

freqCorrect12 = 2*pi*(pulse1.RFamp).^2./(2*centerFreq2);
freqCorrect21 = 2*pi*(pulse2.RFamp).^2./(2*centerFreq1);

phaseCorrect12 = -1/2*cumsum(freqCorrect12)*dt;
phaseCorrect21 = -1/2*cumsum(freqCorrect21)*dt;

pulse1.RFphase = pulse1.RFphase + phaseCorrect21;
pulse2.RFphase = pulse2.RFphase + phaseCorrect12;

% Now, add!
PulseOut = PulseAdd(pulse1, pulse2);

PulseOut.Gx = pulse.Gx;
PulseOut.Gy = pulse.Gy;
PulseOut.Gz = pulse.Gz;