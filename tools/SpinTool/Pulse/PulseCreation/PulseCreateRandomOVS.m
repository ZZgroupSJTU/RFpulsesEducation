function pulse = PulseCreateRandomOVS(passband, totalSW, maxB1, numSteps)
% Creates a random noise pulse for suppression everything outside the
% passband.

numSteps = numSteps + mod(numSteps,2);
numPassSteps = ceil(passband/totalSW*numSteps);
numPassSteps = numPassSteps + mod(numPassSteps, 2);

freqResponse = rand(1,numSteps)-0.5;
numSteps
numPassSteps
passbandIndices = numSteps/2-numPassSteps/2:numSteps/2+numPassSteps/2;
allIndices = 1:numSteps;
stopbandIndices = allIndices; stopbandIndices(passbandIndices)=[];
freqResponse(passbandIndices) = 1;
freqResponse(stopbandIndices) = 0;


figure
dwellTime = 10/numSteps;
sw = 1/dwellTime;
freqAxis = linspace(-sw/2, sw/2, numSteps);
plot(freqAxis, freqResponse)

B1 = ifft((freqResponse));

pulse.tp = 10;
pulse.RFamp = real(B1)./max(abs(B1))*maxB1;
%pulse.RFphase = angle(B1);
pulse.RFphase = zeros(1,numSteps);
pulse.Gx = zeros(1, numSteps);
pulse.Gy = zeros(1, numSteps);
pulse.Gz = zeros(1, numSteps);


