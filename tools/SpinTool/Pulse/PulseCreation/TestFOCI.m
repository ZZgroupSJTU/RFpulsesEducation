clc
close all


numSteps = 1024;
duration = 15;
peakPower = 1;
bandwidth = 18;
truncFactor = 5.2;

normTimeAxis = 2*linspace(0,1, numSteps) - 1;   % Ranges from -1 to 1

flatFraction = 0.3;
numStepsFlat = round(numSteps*flatFraction/2);
numStepsMiddle = numSteps - 2*numStepsFlat;
ampFactor = 6;
t = linspace(-1, 1, numStepsMiddle);
middleShape = ((ampFactor-1))*t.^2 + 1;
shapingFunction = [ampFactor*ones(1, numStepsFlat), middleShape, ampFactor*ones(1, numStepsFlat)];
% fF = 4;
% shapingFunction = cosh(truncFactor*normTimeAxis/2);
% shapingFunction = shapingFunction.*(sech(truncFactor*normTimeAxis/2)>1/fF) + fF.*(sech(truncFactor*normTimeAxis/2)<=1/fF);
% plot(shapingFunction);
% return


freq = -2*pi*bandwidth/(2*truncFactor)*tanh(truncFactor*normTimeAxis);

dwellTime = duration/numSteps;

G = 1;


pulse.tp = duration;
pulse.RFamp = peakPower*sech(truncFactor*normTimeAxis);
pulse.RFphase = cumsum(freq*dwellTime);
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = G*ones(1,numSteps);



pulseFOCI.tp = duration;
RFshape = shapingFunction.*peakPower.*sech(truncFactor*normTimeAxis);
RFshape = RFshape.*(RFshape<=1) + (RFshape>1);
pulseFOCI.RFamp = RFshape;
pulseFOCI.RFphase = cumsum(shapingFunction.*freq*dwellTime);
pulseFOCI.Gx = zeros(1,numSteps);
pulseFOCI.Gy = zeros(1,numSteps);
pulseFOCI.Gz = G*shapingFunction;

pulseFOCI2 = pulseFOCI;
freq2 = -2*pi*bandwidth/(2*truncFactor)*tanh(truncFactor*normTimeAxis) + 25;
pulseFOCI2.RFphase = cumsum(shapingFunction.*freq2*dwellTime);


numSpins = 500;
sampleSize = 500;
xx = linspace(-sampleSize/2, sampleSize/2, numSpins);

offset = 0.7;

pulseFOCIAdded = PulseAddWithDelay(pulseFOCI, pulseFOCI2, 0);

spins0 = InitSpinsRelax(0, numSpins, sampleSize, [0; 0; 1], 1e6, 1e6, 1);
spins0 = ApplyPulseRelax(spins0, pulseFOCIAdded);
spins = InitSpinsRelax(offset, numSpins, sampleSize, [0; 0; 1], 1e6, 1e6, 1);
spins = ApplyPulseRelax(spins, pulseFOCIAdded);
for idx=1:numSpins
    Mz(idx)=spins(idx).M(3);
    Mz0(idx)=spins0(idx).M(3);
end

spins0 = InitSpinsRelax(0, numSpins, sampleSize, [0; 0; 1], 1e6, 1e6, 1);
spins0 = ApplyPulseRelax(spins0, pulseFOCI2);
spins = InitSpinsRelax(offset, numSpins, sampleSize, [0; 0; 1], 1e6, 1e6, 1);
spins = ApplyPulseRelax(spins, pulseFOCI2);
for idx=1:numSpins
    MzFOCI(idx)=spins(idx).M(3);
    MzFOCI0(idx)=spins0(idx).M(3);
end


figure
subplot(3,1,1)
plot(xx,Mz);
hold
plot(xx,Mz0,'r--');
subplot(3,1,2)
plot(xx,MzFOCI);
hold
plot(xx,MzFOCI0,'r--');
subplot(3,1,3)
plot(pulse.RFamp);
hold
plot(pulseFOCI.RFamp, 'r--');