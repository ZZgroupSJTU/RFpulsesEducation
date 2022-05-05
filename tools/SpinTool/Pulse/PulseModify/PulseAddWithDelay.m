function PulseOut = PulseAddWithDelay(pulse1,pulse2,delay)
% SYNTAX: function PulseOut = PulseAdd(pulse1,pulse2,delay);
%
% Adds two RF pulses together, with a certain predefined delay between
% them. The pulses are assumed to have the same number of steps.


% Add zeros before and after
dwellTime = pulse1.tp/numel(pulse1.RFamp);
stepsBefore = round(delay/dwellTime);
stepsAfter = stepsBefore;
if (stepsAfter<0)
    stepsAfter = 0;
end
emptyPulseBefore = PulseCreateZero(delay, stepsBefore);
emptyPulseAfter = PulseCreateZero(delay, stepsAfter);
pulse1 = PulseConcat(emptyPulseBefore, pulse1);
pulse2 = PulseConcat(pulse2, emptyPulseAfter);

PulseOut = PulseAdd(pulse1, pulse2);