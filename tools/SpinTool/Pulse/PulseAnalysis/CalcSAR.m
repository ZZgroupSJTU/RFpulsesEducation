function SAR = CalcSAR(pulse, pulseReference)
% SYNTAX:
%
%    SAR = CalcSAR(pulse, pulseReference)
%
% Calculates the (relative) SAR, in kHz^2*ms, of a pulse.
% The optional input pulseReference is used to normalize the SAR of 
% the input pulse. It can be set to one of two inputs:
%
%    1. Another pulse structure.
%    2. 'rect' - normalize the SAR of the input pulse by that 
%       of a rectangular pulse of equal duration and max amplitude 
%    3. 'ref' - normalize the SAR of the input pulse by the of
%       a 1 ms 0.5 kHz 180 rectangular pulse

numSteps = length(pulse.RFamp);
dwellTime = pulse.tp/numSteps;
SAR = sum(pulse.RFamp.^2*dwellTime);

if (nargin>1)
    if ischar(pulseReference)
        switch lower(pulseReference)
            case 'rect'
                pulseReference = PulseCreateConst(pulse.tp, ...
                                                  numel(pulse.RFamp), ...
                                                  max(pulse.RFamp), ...
                                                  0);
            case 'ref'
                pulseReference = PulseCreateConst(1, 1, 0.5, 0);
            otherwise
                disp('CalcSAR Error: Reference pulse unrecognized. Aborting...');
                SAR = 0;
                beep
                return
        end
    end
    numSteps = length(pulseReference.RFamp);
    dwellTime = pulseReference.tp/numSteps;
    SARReference = sum(pulseReference.RFamp.^2*dwellTime);
    SAR = SAR/SARReference;
end
