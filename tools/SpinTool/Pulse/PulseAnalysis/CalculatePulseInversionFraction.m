function invFrac = CalculatePulseInversionFraction(pulse, minFreq, maxFreq, T1, T2)
% Given an inversion pulse, this function calculates (via simulating
% the Bloch equations) "how much" the pulse inverts within the given
% spectral limits.
%
% Longer explanation: suppose a particular pulse's inversion profile
% between -2 and +2 kHz looks like this:
%
%                             Mz
%                            /|\
%                             |
%              ___            -1           ___
%                 \           |           /
%                  \          |          /
%                   \         |         /
%                    \        |        /
%             -|------\---------------/------|--> Freq. (kHz)
%             -2       \      |      /       2   
%                       \     |     /
%                        \    |    /
%                         \___|___/
%                             |
%                             -(-1)
%                             | 
%
% This obviously imperfect profile of Mz does not invert all of the
% magnetization for three reasons:
% 1. The magnetization in the inversion band does not reach -1
% 2. There is some noninverted band.
% 3. There are transition bands.
% To quantify the "quality of inversion", or the "inversion fraction",
% this function computes the area between the inversion profile and
% the y=1 axis, and divides it by the total area (equal to
% 2*(maxFreq-minFreq)). 
%
%
% Inputs:
%
% Variable Name       Units   Description
% pulse               -       Inversion pulse
% minFreq             kHz     Minimal freq in inversion band
% maxFreq             kHz     Maximal freq in inversion band [1]
% T1                  ms      Longitudinal relaxation
% T2                  ms      Transverse relaxation
% 
% Outputs:
%
% Variable Name       Units   Description
% invFrac             -       The inverted fraction, between 0 and 1.
%                             1 - Perfect inversion profile
%                             0 - Nothing is inverted! Horrible!
% 
% 
% [1] These numbers are arbitrary in a sense that by choosing them to
%     be larger and larger, the inverted fraction becomes smaller and
%     smaller. For a "fair" comparison, one should choose maxFreq-minFreq
%     to equal the theoretical pulse bandwidth (e.g. ~ 1/duration for
%     simple Fourier pulses).

numPoints = 300;
initMag = [0; 0; 1];
linPhaseCorrect = 0;
wrapPhase = 0;

[~ , ~ , Mz, ~, ~, ~] = CalcPulseFreqResponseWithRelax(...
                                 pulse, ...
                                 T1, ...
                                 T2, ...
                                 minFreq, ...
                                 maxFreq, ...
                                 numPoints, ...
                                 initMag, ...
                                 linPhaseCorrect, ...
                                 wrapPhase);

% If Mz = -1 everywhere, sum(Mz) = -numPoints.
% If Mz = +1 everywhere, sum(Mz) = numPoints
invFrac = (numPoints-sum(Mz))./(2*numPoints);   