%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% This script creates a movie demonstrating the concept of spin echoes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clc
clear all
close all

cs = randn(1,2500)*0.1;
numSpinsPerCS = 1;
initialMag = [0;0;1];
initialMag = initialMag./norm(initialMag);
sampleSize = 1;
T1 = 1e9;
T2 = 1e9;
spins = InitSpinsRelax (cs, numSpinsPerCS, sampleSize, initialMag, T1, T2);

delayTime = 20; % ms
delayNumSteps = 200;
numRFSteps = 20;
RFDuration = 0.1; % ms
seq = {{'rect', 90, 0, RFDuration, numRFSteps}, {'delay', delayTime, delayNumSteps}, {'rect', 180, 90, RFDuration, 150}, {'delay', delayTime*2, delayNumSteps*2}};

numFrames = 250;
numPlottedSpins = 15;
maxMovieDuration = CalcSeqTotalTime(seq);
isCloseWindow = 0;
isInvertedColor = 1;
[mov1, spins, fid1] = MovieBlochIsochromatsSeq(...
    seq, ...
    'spins', spins, ...
    'numFrames', numFrames, ...
    'numPlottedSpins', numPlottedSpins, ...
    'isPlotRF', false, ...
    'isCloseWindow', false, ...
    'isPlotFID', true, ...
    'isInvertedColor', true, ...
    'maxMovieDuration', maxMovieDuration);

