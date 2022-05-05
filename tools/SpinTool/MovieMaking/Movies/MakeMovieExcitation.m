% ========================================================================
% Create a movie showing the trajectory of a spin during an excitation
% pulse. Amongst other things, I want to see how the spin acquires the
% linear phase term exp(-i*w*T/2) as it precesses.
% ========================================================================

clc
close all
clear all

flipAngle = 90; % deg
pulseDuration = 3.84; % ms
[pulse, pulseHdr] = PulseReadSiemensInclude('SLR90.h', pulseDuration, flipAngle);
% pulse.RFphase = pulse.RFphase - pi/2;
seq = {pulse};

% Create the spin
csVec = 0.15;
numSpins = 1;
sampleSize = 1;
initMag = [0;0;1];
T1 = 1e6;
T2 = 1e6;
eqMag = 1;
spins = InitSpinsRelax(csVec, numSpins, sampleSize, initMag, T1, T2, eqMag);

% [movieSpins, movieTimeAxis, fid]=MovieBlochIsochromatsSeq(seq,spins,numFrames,numPlottedSpins, plotRF, isCloseWindow, isInvertedColor, isPlotFID, isPlotSpinTrajectories, maxMovieDuration, spinColorVec)
numFrames = 300;
numPlottedSpins = 1;
plotRF = 0;
isCloseWindow = 1;
isInvertedColor = 0;
isPlotFID = 0;
isPlotSpinTrajectories = 1;
maxMovieDuration = 4;
spinColorVec{1} = [0 0 1];
[movieSpins, movieTimeAxis, fid] = MovieBlochIsochromatsSeq(seq, spins, numFrames, numPlottedSpins, plotRF, isCloseWindow, isInvertedColor, isPlotFID, isPlotSpinTrajectories, maxMovieDuration, spinColorVec);
