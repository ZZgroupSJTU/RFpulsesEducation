function pulse = PulseCreateChirp(initialFreq, finalFreq, pulseDuration, numSteps, sampleLength, whichAxis)
% Description: The current function creates a pulse structure for 
% a chirped, linearly swept, homogeneous excitation pulse. 
%
% Inputs:
%
% Variable Name   Units   Description
% initialFreq     kHz     Initial excitation frequency
% finalFreq       kHz     Final excitation frequency
% pulseDuration   ms      Pulse duration
% numSteps        -       Number of excitation points
% sampleLength    mm      Sample length (used to set gradient). Optional
% whichAxis       char    Equal to 'x', 'y', 'z'. Optional
%
% Outputs:
%
% Variable Name   Units   Description
% pulse           -       OupulseDurationut pulse structure.

% Create time vector
excTimeVec = linspace(0,pulseDuration,numSteps);

% Calculate the sweep rate (often denoted R in many papers)
bandwidth = finalFreq - initialFreq;
sweepRate = bandwidth/pulseDuration;

% Set pulse parameters. 
pulse.tp = pulseDuration;
pulse.RFamp = wurst(40,numSteps); % Use a wurst-shaped RF amplitude to smooth the response of the pulse. 
pulse.RFphase = 2*pi*(initialFreq*excTimeVec + sweepRate*excTimeVec.^2 / 2); % Quadratic phase profile
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);

% Set chirp's amplitude for a 90-deg. excitation. For details, see
% appendix in Shrot & Frydman, JMR 172:179-190 (2005).
pulse.RFamp = pulse.RFamp./max(pulse.RFamp)*0.267*sqrt(abs(sweepRate));

% Set an appropriate field gradient
if nargin>4
    gradientStrength = bandwidth/sampleLength;   % kHz/mm
    switch whichAxis
        case 'x',
            pulse.Gx = gradientStrength.*ones(1,numSteps);
        case 'y',
            pulse.Gy = gradientStrength.*ones(1,numSteps);
        case 'z',
            pulse.Gz = gradientStrength.*ones(1,numSteps);
        otherwise,
            disp('Error in PulseCreateChirp: axis undefined. Aborting!');
            beep
            return
    end;
end

