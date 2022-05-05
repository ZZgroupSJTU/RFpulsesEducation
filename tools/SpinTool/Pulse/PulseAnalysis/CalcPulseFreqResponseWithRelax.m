function [Mx, My, Mz, phaseMxy, magMxy, freqAxis] = CalcPulseFreqResponseWithRelax(...
                                                     pulse, ...
                                                     T1, ...
                                                     T2, ...
                                                     freqMin, ...
                                                     freqMax, ...
                                                     numResPoints, ...
                                                     initMag, ...
                                                     linPhaseCorrect, ...
                                                     wrapPhase)
% SYNTAX: PulseFreqResponse(pulse,freqMin,freqMax,N,initMag,linPhaseCorrect,wrapPhase);
%
% The current function plots the frequency response of the input pulse from
% frequency freqMin to frequency freqMax, with N steps inbetween. It assumes
% all the spins start from the initial condition initMag which is supplied (and
% is assumed to be a 3x1 vector).
%
% Note that for such a plot, the gradient fields in the pulse object aren't
% used, only pulse.RFamp, pulse.RFphase and pulse.tp play a role.
%
% Optional parameters:
%     linPhaseCorrect - linear phase correction
%     wrapPhase   - set to 1 for phase wrapping (interval [0,2*pi])

if nargin<6
    linPhaseCorrect = 0;
    wrapPhase = 0;
end;

if nargin<7
    wrapPhase = 0;
end;

% Step I: Create a spin structure
freqAxis = linspace(freqMin, freqMax, numResPoints);  % kHz
spins = InitSpinsRelax(freqAxis, 1, 0, initMag, T1, T2, 1);

% Apply pulse
spins = ApplyPulseRelax(spins,pulse);

% Extract output magnetization
M = zeros(3, numResPoints);  % Pre-allocate
for k=1:numResPoints
    M(:,k) = spins(k).M;
end;

% Calculate magMxy, phaseMxy, Mx, My (with linPhaseCorrect taken into account)
magMxy = sqrt(M(1,:).^2 + M(2,:).^2);
phaseMxy = phase(M(1,:)+1i*M(2,:)) + freqAxis*linPhaseCorrect;
if wrapPhase==1
    phaseMxy = mod(phaseMxy,2*pi);
end;
Mx = M(1,:); 
My = M(2,:); 
Mz = M(3,:); 
Mxy = Mx + 1i*My;