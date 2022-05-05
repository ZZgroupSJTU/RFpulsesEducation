% ========================================================================
%
% This script creates the movies for the talk given in ELSC in Aviv 
% Mezer's group.
% 
% #    Demonstration Purpose       Description
% 1    Larmor precession           Single spin excite-acquire
% 2    T2* decay                   Pulse-acquire on distribution of spins
%                                  showing the signal as it decays
% 3    Spin-Echo                   Pulse-delay-180-delay on distribution 
%                                  of spins showing the signal as it first
%                                  decays (T2*) and then builds up
% 4    Multiple echoes             Echoes arising from multiple 180s
%                                  strung one after the other
% 5    Echoes from 3 pulses        pulse-delay-pulse-delay-pulse-delay
%                                  showing all 5 echos + 3 FIDs forming
%                                  from 3 general pulses (I'll take 
%                                  110-108-108)
% ========================================================================

clc
clear all
close all

numSpinsPerCS = 1;
sampleSize = 1;
T1 = 1e9;
T2 = 1e9;
initialMag = [0;0;1];
initialMag = initialMag./norm(initialMag);

%% ======================================================================= 
% 1. Larmor Precession, single spin
% ======================================================================== 

cs = 0.25;
spins = InitSpinsRelax (cs, numSpinsPerCS, sampleSize, initialMag, T1, T2);

delayTime = 20; % ms
delayNumSteps = 200;
numRFSteps = 20;
RFDuration = 0.1; % ms
seq = {{'rect', 90, 0, RFDuration, numRFSteps}, {'delay', delayTime, delayNumSteps}};

[mov1, spins1, fid1] = MovieBlochIsochromatsSeq(...
    seq, ...
    'spins', spins, ...
    'numFrames', 200, ...
    'numPlottedSpins', 1, ...
    'isPlotRF', 1, ...
    'isCloseWindow', 0, ...
    'isPlotFID', true, ...
    'isInvertedColor', 0, ...
    'isPlotSpinsXY', false, ...
    'maxMovieDuration', CalcSeqTotalTime(seq));


%% ======================================================================= 
% 2. T2* decay
% ======================================================================== 

cs = randn(1,2500)*0.1;
histogram(cs);
spins = InitSpinsRelax (cs, numSpinsPerCS, sampleSize, initialMag, T1, T2);

delayTime = 20; % ms
delayNumSteps = 200;
numRFSteps = 20;
RFDuration = 0.1; % ms
seq = {{'rect', 90, 0, RFDuration, numRFSteps}, {'delay', delayTime, delayNumSteps}};

[mov2, spins2, fid2] = MovieBlochIsochromatsSeq(...
    seq, ...
    'spins', spins, ...
    'numFrames', 200, ...
    'numPlottedSpins', 25, ...
    'isPlotRF', 1, ...
    'isCloseWindow', 0, ...
    'isPlotFID', 1, ...
    'isInvertedColor', 0, ...
    'isPlotSpinsXY', true, ...
    'maxMovieDuration', CalcSeqTotalTime(seq));


%% ======================================================================= 
% 3. Spin echoes (two of them)
% ======================================================================== 

cs = randn(1,2500)*0.1;
spins = InitSpinsRelax (cs, numSpinsPerCS, sampleSize, initialMag, T1, T2);

delayTime = 20; % ms
delayNumSteps = 200;
numRFSteps = 20;
RFDuration = 0.01; % ms
seq = {{'rect', 90, 0, RFDuration, numRFSteps}, {'delay', delayTime, delayNumSteps}, ...
       {'rect', 180, 90, RFDuration, 150}, {'delay', delayTime*2, delayNumSteps*2}, ...
       {'rect', 180, 90, RFDuration, 150}, {'delay', delayTime*2, delayNumSteps*2}};

[mov3, spins3, fid3] = MovieBlochIsochromatsSeq(...
    seq, ...
    'spins', spins, ...
    'numFrames', 250, ...
    'numPlottedSpins', 25, ...
    'isPlotRF', 1, ...
    'isCloseWindow', 0, ...
    'isPlotFID', 1, ...
    'isInvertedColor', 1, ...
    'isPlotSpinsXY', true, ...
    'maxMovieDuration', CalcSeqTotalTime(seq));

% vidObj=VideoWriter('03 - Spin Echoes.avi'); 
% vidObj.FrameRate=10; 
% open(vidObj); 
% vidObj.writeVideo(mov3);
% close(vidObj);

%% ======================================================================= 
% 3b. Spin echoes (just one)
% ======================================================================== 

cs = randn(1,1000)*0.1;
spins = InitSpinsRelax (cs, numSpinsPerCS, sampleSize, initialMag, T1, T2);

delayTime = 20; % ms
delayNumSteps = 200;
numRFSteps = 20;
RFDuration = 0.01; % ms
seq = {{'rect', 90, 0, RFDuration, numRFSteps}, {'delay', delayTime, delayNumSteps}, ...
       {'rect', 180, 90, RFDuration, 150}, {'delay', delayTime*2, delayNumSteps*2}};

[mov3b, spins3b, fid3b] = MovieBlochIsochromatsSeq(...
    seq, ...
    'spins', spins, ...
    'numFrames', 250, ...
    'numPlottedSpins', 15, ...
    'isPlotRF', 1, ...
    'isCloseWindow', 0, ...
    'isPlotFID', 1, ...
    'isInvertedColor', false, ...
    'isPlotSpinsXY', true, ...
    'maxMovieDuration', CalcSeqTotalTime(seq));

vidObj=VideoWriter('03b - Spin Echo.avi'); 
vidObj.FrameRate=10; 
open(vidObj); 
vidObj.writeVideo(mov3b);
close(vidObj);

%% ======================================================================= 
% 4. Reduced flip angle spin echo
% ======================================================================== 

cs = randn(1,2500)*0.1;
spins = InitSpinsRelax (cs, numSpinsPerCS, sampleSize, initialMag, T1, T2);

delayTime = 20; % ms
delayNumSteps = 200;
numRFSteps = 20;
RFDuration = 0.01; % ms
seq = {{'rect', 90, 0, RFDuration, numRFSteps}, {'delay', delayTime, delayNumSteps}, ...
       {'rect', 90, 90, RFDuration, 150}, {'delay', delayTime*2, delayNumSteps*2}};

[mov4, spins4, fid4] = MovieBlochIsochromatsSeq(...
    seq, ...
    'spins', spins, ...
    'numFrames', 250, ...
    'numPlottedSpins', 25, ...
    'isPlotRF', 1, ...
    'isCloseWindow', 0, ...
    'isPlotFID', 1, ...
    'isInvertedColor', 1, ...
    'maxMovieDuration', CalcSeqTotalTime(seq));


%% ======================================================================= 
% 5. Three-pulse echoes
% ======================================================================== 

cs = randn(1,2500)*0.1;
spins = InitSpinsRelax (cs, numSpinsPerCS, sampleSize, initialMag, T1, T2);

D1 = 7; % ms
D2 = 23;
D3 = 64;
delayNumSteps = 200;
numRFSteps = 20;
RFDuration = 0.01; % ms
seq = {{'rect', 110, 0, RFDuration, numRFSteps}, {'delay', D1, delayNumSteps}, ...
       {'rect', 108, 90, RFDuration, 150}, {'delay', D2, delayNumSteps*2}, ...
       {'rect', 108, 90, RFDuration, 150}, {'delay', D3, delayNumSteps*2}};

[mov5, spins5, fid5] = MovieBlochIsochromatsSeq(...
    seq, ...
    'spins', spins, ...
    'numFrames', 250, ...
    'numPlottedSpins', 25, ...
    'isPlotRF', 1, ...
    'isCloseWindow', 0, ...
    'isPlotFID', 1, ...
    'isInvertedColor', 1, ...
    'maxMovieDuration', CalcSeqTotalTime(seq));