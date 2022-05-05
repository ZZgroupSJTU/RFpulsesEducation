function B1_RMS = CalcPulseRMSB1(pulse)
% Calculates the RMS of the RF amplitude (in kHz). 

% B1_RMS = std(pulse.RFamp);
B1_RMS = sqrt(sum(pulse.RFamp.^2)/numel(pulse.RFamp));

