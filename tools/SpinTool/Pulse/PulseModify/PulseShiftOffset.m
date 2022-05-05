function pulse = PulseShiftOffset(pulse, offset, offsetAxis)
% SYNTAX: pulse = PulseShiftOffset(pulse, offset, offsetType)
%
% Adds a linear offset to the given pulse.
%
% Input Variables
% Variable Name       Units            Description
% pulse 
% offset              See Below        Amount of shift
% offsetAxis          'x', 'y', 'z'    Optional argument. Can be either
%                                      'khz' or 'mm', to indicate whether
%                                      the offset should be in kHz or mm.
%                                      If omitted, defaults to 'khz'. If
%                                      set to 'mm', the 'offsetAxis'
%                                      variable must be specified as well
%                                      to indicate
%   

if (nargin<3)
    offsetKHz = offset;
end
if (nargin==3)
    switch lower(offsetAxis)
        case 'z'
            [maxGrad, maxIdx] = max(abs(pulse.Gz));
            offsetGrad = maxGrad*sign(pulse.Gz(maxIdx));
        case 'y'
            [maxGrad, maxIdx] = max(abs(pulse.Gy));
            offsetGrad = maxGrad*sign(pulse.Gz(maxIdx));
        case 'x'
            [maxGrad, maxIdx] = max(abs(pulse.Gx));
            offsetGrad = maxGrad*sign(pulse.Gz(maxIdx));
        otherwise
            disp('Warning: offsetAxis not equal to x, y or z in PulseShiftOffset! Setting offsetGrad to 0.');
            offsetGrad = 0;
    end
    if (offsetGrad==0)
        disp('Warning: offset gradient zero in PulseShiftOffset. Setting offset to 0!');
        offsetKHz = 0;
    else
        offsetKHz = offsetGrad*offset; % kHz/mm * mm = kHz
    end
end


numSteps = length(pulse.RFamp);
dt = pulse.tp/numSteps;
timeAxis = [0:dt:dt*(numSteps-1)];

%               Original phase    Linear phase            Ensures phase=0 at center of pulse
pulse.RFphase = pulse.RFphase - 2*pi*offsetKHz*timeAxis + pi*offsetKHz*pulse.tp;