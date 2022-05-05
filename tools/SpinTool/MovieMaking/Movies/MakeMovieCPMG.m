% ========================================================================
% Creates several movies illustrating the dynamics of spins under an
% external magnetic field (no relaxation).
% ========================================================================

clear all
close all
clc

% ========================================================================
% Precession around a constant magnetic field
% ========================================================================

pulseDuration = 30; % ms
spinPos = [0;0;0]; % mm
numFramesToSkip = 1;
addAnnotations = 0;
isNormGlobalB = 0.9;
isShowNormalPlane = 0;
isPlotSpinTrajectory = 1;
numPulseSteps = 200;

csVec = randn(1,491)*0.3;
numSpinsPerShift = 1;
sampleSize = 10;
initialMag = [0;0;1];
T1 = 1e6;
T2 = 1e6;
eqMag = 1;
spins = InitSpinsRelax(csVec, numSpinsPerShift, sampleSize, initialMag, T1, T2, eqMag);

B1Inho = 1; % Set to 1 for no inhomogeneity
Ph90 = 0; 
Ph180 = Ph90 + 90;
T90 = 0.01; % ms
T180 = 0.01; % ms
delayNumSteps = 1001;
delayTime = 3;
% seqCPMG = {{'rect', 90*B1Inho, Ph90, T90, 150}, {'delay', delayTime, delayNumSteps}, {'rect', 180*B1Inho, Ph180, T180, 150}, {'delay', delayTime*2, delayNumSteps*2},  {'rect', 180*B1Inho, Ph180, T180, 150}, {'delay', delayTime*2, delayNumSteps*2}, {'rect', 180*B1Inho, Ph180, T180, 150}, {'delay', delayTime*2, delayNumSteps*2}, {'rect', 180*B1Inho, Ph180, T180, 150}, {'delay', delayTime, delayNumSteps}};
seqCPMG = {{'rect', 90*B1Inho, Ph90, T90, 150}, {'delay', delayTime, delayNumSteps}, {'rect', 180*B1Inho, Ph180, T180, 150}, {'delay', delayTime*2, delayNumSteps*2}};
% seqCPMG = {{'rect', 90*B1Inho, Ph90, T90, 150}, {'delay', delayTime, delayNumSteps}, {'rect', 180*B1Inho, Ph180, T180, 150}, {'delay', delayTime*2, delayNumSteps*2},  {'rect', 180*B1Inho, Ph180, T180, 150}, {'delay', delayTime*2, delayNumSteps*2}};

maxMovieDuration = CalcSeqTotalTime(seqCPMG);
isInvertedColors = 0;
%camPosition = [
%    PlotBloch
plotRF = 1;
numPlottedSpins = 15;
isCloseWindow = 0;
numFrames = 400;
isPlotFID = 1;
[outputMovieCPMG, movieTimeAxis, fid] = MovieBlochIsochromatsSeq (seqCPMG, spins, numFrames, numPlottedSpins, plotRF, isCloseWindow, isInvertedColors, isPlotFID, maxMovieDuration);
% figure
% plot(movieTimeAxis, imag(fid)/numPlottedSpins)