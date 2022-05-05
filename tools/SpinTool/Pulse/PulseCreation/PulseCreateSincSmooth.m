function pulse = PulseCreateSincSmooth(numLobes, duration, numSteps, flipAngle, wurstParameter)

pulse = PulseCreateSinc(numLobes, duration, numSteps, flipAngle);
pulse.RFamp = pulse.RFamp.*wurst(wurstParameter, numSteps);
