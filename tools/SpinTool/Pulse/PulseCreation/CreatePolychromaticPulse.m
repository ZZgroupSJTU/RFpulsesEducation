function pulseTrain = CreatePolychromaticPulse(pulse, centerFreqs, pulseDelay, phaseVector, ampVector)
% Description: uses a given pulse profile to generate a 
% polychromatic pulse, made up of copies of the given 
% pulse spaced a certain delay apart, with given center 
% frequencies.
%
% Inputs:
%
% Variable Name       Units   Description
% pulse               -       Given pulse
% centerFreqs         kHz     Vector of center frequencies 
% pulseDelay          ms      Delay between pulse copies
% phaseVector         rad     Relative phases of pulses,
%                             (e.g. [pi 0 pi 0] for 4 pulses)
% ampVector           a.u.    Relative amplitude scaling of pulses
%                             (e.g. [1 1 1 1] for 4 pulses means nothing)
%
% Outputs:
%
% Variable Name       Units   Description
% pulseTrain          -       Pulse train


% ********************************
% Retrieve single pulse parameters
% ********************************

pulseDuration = pulse.tp;
pulseNumSteps = length(pulse.RFamp);
dwellTime = pulseDuration/pulseNumSteps;

% Make sure the delay time is a multiple of the basic dwell time
pulseDelay = ceil(pulseDelay/dwellTime)*dwellTime;

% ******************************
% Compute pulse train properties
% ******************************

% Number of pulses in train
numPulses = length(centerFreqs);

% Total pulse duration 
trainDuration = (numPulses-1)*pulseDelay + pulseDuration;

% Total number of steps
trainNumSteps = ceil(trainDuration/dwellTime);

% Initialize final pulse
pulseTrain = PulseCreateZero(trainDuration, trainNumSteps);

% Add 'em up!
for curPulseIdx = 1:numPulses
    % Frequency shift the sech pulse to the correct frequency center
    shiftedPulse = PulseShiftOffset(pulse, centerFreqs(curPulseIdx));
    % Adjust global phase and amplitude. Wrap phase to [0 2pi]
    shiftedPulse.RFamp = shiftedPulse.RFamp*ampVector(curPulseIdx);
    shiftedPulse.RFphase = shiftedPulse.RFphase + phaseVector(curPulseIdx);
    shiftedPulse.RFphase = mod(shiftedPulse.RFphase, 2*pi);
    % First, create a pulse padded with 0s before and after to reflect the
    % delays. Make sure the delay divides by the dwell time; if not,
    % add more pulse steps.
    timeBefore = (curPulseIdx-1)*pulseDelay;
    stepsBefore = round(timeBefore/dwellTime);
    timeAfter = (numPulses - curPulseIdx)*pulseDelay;
    stepsAfter = trainNumSteps -  stepsBefore - pulseNumSteps;
    if (stepsAfter<0)
        stepsAfter = 0;
    end
    emptyPulseBefore = PulseCreateZero(timeBefore, stepsBefore);
    emptyPulseAfter = PulseCreateZero(timeAfter, stepsAfter);
    % Next, create a frequency shifted version of the 
    curPulse = PulseConcat(emptyPulseBefore, shiftedPulse);
    curPulse = PulseConcat(curPulse, emptyPulseAfter);
    % Now add the current pulse to the other ones
    pulseTrain = PulseAdd(pulseTrain, curPulse);
end