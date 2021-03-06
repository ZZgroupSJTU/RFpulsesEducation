function pulseout = PulseCalibrate90(pulse,offset,calibrationParameter)
% SYNTAX: function pulseout = PulseCalibrate90(pulse,offset,num_of_iterations)
%
% The current function calibrates the power of the input pulse such that the
% resulting spin at the offset (in kHz) will have reached the xy plane, i.e. 
% Mz = 0. It assumes the offset starts from the state Mx = My = 0, Mz = 1.
%
% calibrationParameter is a string, set to either 'amplitude' or 'duration'
% indicating what exactly should be calibrated. Is not entered, "amplitude"
% is assumed.
%
% num_of_iterations is the number of iterations the program carries out before
% exiting. More iterations = more time = more accuracy.
%
% The function returns the same pulse as the input, only with calibrated power.


pulseout = pulse;

switch (calibrationParameter)
    case 'duration' % Calibrate pulse's duration
        initialDurationGuess = 0.01;
        duration = fzero(@(x) Calc_Mz_Duration(x,pulse,offset),initialDurationGuess);
        pulseout.tp = duration;
    otherwise % Calibrate pulse's power
        initialRFGuess = 0.01;
        RFpower = fzero(@(x) Calc_Mz_Power(x,pulse,offset),initialRFGuess);
        pulseout.RFamp = pulseout.RFamp./max(abs(pulseout.RFamp)).*RFpower;
end





% ------------------------ ADDITIONAL FUNCTIONS -----------------------



function Mz=Calc_Mz_Power(RFpower, pulse,offset)

if (RFpower<0)
    Mz = 1;
else
    % Adjust RFpower
    pulse.RFamp = pulse.RFamp./max(abs(pulse.RFamp)).*RFpower;
    % Create input spin
    spin = InitSpinsRelax(offset, 1, 1, [0; 0; 1], 1e6, 1e6, 1);
    % Simulate
    SpinOut = ApplyPulseRelax(spin,pulse);
    Mz = SpinOut(1).M(3);
end

function Mz=Calc_Mz_Duration(duration, pulse,offset)

if (duration<0)
    Mz = 1;
else
    % Adjust RF duration (in ms)
    pulse.tp = duration;
    % Create input spin
    spin = InitSpinsRelax(offset, 1, 1, [0; 0; 1], 1e6, 1e6, 1);
    % Simulate
    SpinOut = ApplyPulseRelax(spin,pulse);
    % Save result
    Mz = SpinOut(1).M(3);
end