close all
clc

maxPeakPower = 1; % kHz
powerRatio = 1;
hadamardOrder = 1;
FOV = 1; % mm
fieldInhomogeneity = 0;
fillingFactor = 1;
GMR = 42; % kHz/mT
sampleSize = 10; % mm
maxGrad = 40; % mT/m
truncFactor = 5.2;
adiabaticityFactor = 4;
minDwellTime = 0.001;

% *******************************
% Compute single pulse parameters
% *******************************

% Pulse power
maxPeakPowerSinglePulse = maxPeakPower/powerRatio;

% Voxel size (in meters)
numVoxels = hadamardOrder;
voxelSize = FOV/hadamardOrder;

% Duration and bandwidth. 
% First, compute constant of proportionality, C, between both: T = C*G
pulseGrad = 2*(maxPeakPowerSinglePulse-fieldInhomogeneity)/(fillingFactor*GMR*voxelSize);

if (pulseGrad>maxGrad)
    disp(['Maximal gradient value exceeded (',num2str(pulseGrad),' mT/m), resetting to maximal value (',num2str(maxGrad),' mT/m).']);
    pulseGrad = maxGrad;
end

% Find the corresponding pulse duration, given the adiabaticity factor
pulseDuration = 4*truncFactor*adiabaticityFactor/(2*pi*fillingFactor*GMR*voxelSize*pulseGrad);

% Bandwidth per voxel (kHz)
voxelBW = GMR*pulseGrad*voxelSize;

% FOV Bandwidth (kHz) = num of Voxels * voxelBW
FOVBW = pulseGrad*GMR*FOV; % kHz

% Pulse bandwidth
pulseBW = voxelBW*fillingFactor;

% Sample bandwidth and dwell time (ms)
sampleBW = GMR*pulseGrad*sampleSize;
dwellTime = 1/sampleBW;  

% Compute number of pulse steps. Round up, if necessary, and recalculate
% dwell time accordingly
numSteps = ceil(pulseDuration/dwellTime);
dwellTime = pulseDuration/numSteps; 

% Verify dwell time isn't too short
if (dwellTime<minDwellTime)
    disp(['Dwell time = ',num2str(dwellTime*1000),' us is too short - aborting.']);
    return
end

% ********************
% Compute single pulse 
% ********************

pulseSech = PulseCreateAdiabatic(...
          pulseBW, ...
          truncFactor, ...
          maxPeakPowerSinglePulse, ...
          pulseDuration, ...
          numSteps, ...
          'sech');

pulseSech8 = PulseCreateAdiabatic(...
          pulseBW, ...
          truncFactor, ...
          maxPeakPowerSinglePulse*0.5, ...
          pulseDuration, ...
          numSteps, ...
          'sech8');

      
% ****************
% Simulate results
% ****************

T1 = 10000;
T2 = 10000;
freqMin = -sampleSize/2*pulseGrad*GMR*0.3;
freqMax = sampleSize/2*pulseGrad*GMR*0.3;
numPoints = 500;
initMag = [0; 0; 1];

[~, ~, MzSech, ~, ~, ~] = CalcPulseFreqResponseWithRelax(pulseSech, T1, T2, freqMin, freqMax, numPoints, initMag, 0, 0);
[~, ~, MzSech8, ~, ~, freqAxis] = CalcPulseFreqResponseWithRelax(pulseSech8, T1, T2, freqMin, freqMax, numPoints, initMag, 0, 0);

figure
plot(pulseSech.RFamp);
hold
plot(pulseSech8.RFamp,'r');

figure
plot(freqAxis,  MzSech);
hold
plot(freqAxis,  MzSech8, 'r');
hold off

