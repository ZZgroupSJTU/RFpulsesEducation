clc
close all

%excitationBW = 2;
FOV = 1500;
peakB1 = 0.8;
minSlice = 10; % mm
isExport = 1;
[pulse, BW] =PulseCreateSLR90(FOV, peakB1, minSlice, 'linearr12', isExport);




T1 = 1e6;
T2 = 1e6;
numResPoints = 1000;
initMag = [0; 0; 1];
linPhaseCorrect = 0;
wrapPhase = 0;
FOVBW = FOV/minSlice*4;
[Mx, My, Mz, phaseMxy, magMxy, freqAxis] = CalcPulseFreqResponseWithRelax(...
                                                     pulse, ...
                                                     T1, ...
                                                     T2, ...
                                                     -FOVBW, ...
                                                     FOVBW, ...
                                                     numResPoints, ...
                                                     initMag, ...
                                                     linPhaseCorrect, ...
                                                     wrapPhase);
                                                
plot(magMxy);                                                 
gradFor10mm = BW/(GetGyromagneticRatio('1h')*0.01);
pulse.tp;
sum(magMxy(1:425))*2/sum(magMxy);
%ExportPulseToSiemens(pulse, 'SLR90', 'SLR90', BW, pi/2, ['SLR Linear 90. Fractional |Mxy| area outside BW ~ 0.006. Peak = 1 kHz @ 3.2 ms. FOV = ', num2str(FOVBW/BW)]);
%ExportPulseToSiemens(pulse, 'SLR90', 'SLR90', BW, pi/2, ['SLR Linear 90. Fractional |Mxy| area outside BW ~ 0.6\%']);

