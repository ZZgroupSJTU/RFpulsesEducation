function [pulseTrain, pulseInfo] = CreatePolychromaticExcitationPulse(basePulse, numBands, fillingFactor, maxB1)
% Creates a polychromatic pulse excitating multiple, equi-spaced lobes,
% by using the input basePulse as a reference. The routine automatically
% calibrates the required parameters.


% ------------------------------------------------------------------------
% Calibrate pulse's duration for 90
% ------------------------------------------------------------------------

basePulse.RFamp = basePulse.RFamp./max(abs(basePulse.RFamp))*maxB1/numBands;
calibrationOffset = 0; % Which frequency to use to determine 90
basePulse = PulseCalibrate90(basePulse,calibrationOffset,'duration');


% ------------------------------------------------------------------------
% Calculate pulse's excitable bandwidth
% ------------------------------------------------------------------------

% Calculate pulse's frequency response. It is safe to assume that 
% the pulse's frequency bandwidth is not larger than 20/(it's duration).
% If we want 100 points of resolution for the pulse's excitatory shape,
% we then need to take 20 times take to account for the larger simulated
% bandwidth, and hence numResPoints = 2000.
T1 = 1e6;
T2 = 1e6;
numResPoints = 2000;
initMag = [0; 0; 1];
[~, ~, Mz, ~, ~, freqAxis] = CalcPulseFreqResponseWithRelax(basePulse, T1, T2, -10/basePulse.tp, 10/basePulse.tp, numResPoints, initMag, 0, 0);

% Calculate width of excitation lobe @ Mz = 0.05 
indicesMz = find(Mz<=0.05);
bandwidthAtBase = freqAxis(indicesMz(end)) - freqAxis(indicesMz(1));

% To make sure the excitation profile is "good", we also check for the FWHM
% and compare it to the bandwidthAtBase
indicesMz = find(Mz<=0.5);
FWHM = freqAxis(indicesMz(end)) - freqAxis(indicesMz(1));

if (FWHM>2*bandwidthAtBase)
    disp('Caution in CreatePolychromaticExcitationPulse: basePulse might be malformed? FWHM seems too wide!');
end

% ------------------------------------------------------------------------
% Form the polychromatic pulse
% ------------------------------------------------------------------------

pulseDelay = 0;

% Determine center frequencies of excitation lobes
distanceBetweenCenters = bandwidthAtBase/fillingFactor;
centerFrequencies = linspace(-(numBands-1)*distanceBetweenCenters/2, (numBands-1)*distanceBetweenCenters/2, numBands);

pulseTrain = CreatePolychromaticPulse(basePulse, centerFrequencies, pulseDelay);


totalBW = distanceBetweenCenters*numBands;
[~, ~, MzTotal, ~, MxyTotal, freqAxisTotal] = CalcPulseFreqResponseWithRelax(pulseTrain, T1, T2, -totalBW/2*1.5, totalBW/2*1.5, numResPoints, initMag, 0, 0);

% ------------------------------------------------------------------------
% Prepare pulse information
% ------------------------------------------------------------------------

pulseInfo.totalBW = totalBW;
pulseInfo.FWHM = FWHM;
pulseInfo.Mz = Mz;
pulseInfo.bandwidthPerBand = distanceBetweenCenters;
pulseInfo.freqAxis = freqAxis;
pulseInfo.MxyTotal = MxyTotal;
pulseInfo.MzTotal = MzTotal;
pulseInfo.freqAxisTotal = freqAxisTotal;