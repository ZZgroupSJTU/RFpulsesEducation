function pulse = PulseCreateSechALEPH(peakPower, adiabaticity, numSteps)

eta = 0.8;
kappa = 1.7;
truncFactor = 5.3;

B1Threshold = peakPower*(1-adiabaticity)

AFPRc = (B1Threshold/eta)^2;
AFPBW = kappa*sqrt(AFPRc)
AFPDuration = truncFactor*AFPBW/AFPRc;


pulse = PulseCreateSech(AFPBW, truncFactor, peakPower, AFPDuration, numSteps, 0);