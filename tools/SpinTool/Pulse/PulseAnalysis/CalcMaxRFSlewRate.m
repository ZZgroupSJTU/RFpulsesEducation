function maxRFSlewRate = CalcMaxRFSlewRate(pulse)

numSteps = numel(pulse.RFamp);
dt = pulse.tp/numSteps;
RF = pulse.RFamp.*exp(1i*pulse.RFphase);
maxRFSlewRate = max(abs(diff(RF))*dt);  % kHz/ms


