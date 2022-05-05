function PulseOut = PulseAdd(pulse1,pulse2,pulseSpatial)
% SYNTAX: function PulseOut = PulseAdd(pulse1,pulse2[,pulseSpatial]);
%
% Adds two RF pulses together. Pulses are assumed to have the same duration
% and number of steps. Gradients are set to 0.
%

fieldVec1x = pulse1.RFamp.*cos(pulse1.RFphase);
fieldVec1y = pulse1.RFamp.*sin(pulse1.RFphase);
fieldVec2x = pulse2.RFamp.*cos(pulse2.RFphase);
fieldVec2y = pulse2.RFamp.*sin(pulse2.RFphase);
fieldx = fieldVec1x + fieldVec2x;
fieldy = fieldVec1y + fieldVec2y;
fieldAmp = sqrt(fieldx.^2 + fieldy.^2);
fieldPhase = atan2(fieldy, fieldx);
numSteps = length(pulse1.RFamp);


PulseOut.tp = pulse1.tp;
PulseOut.RFamp = fieldAmp;
PulseOut.RFphase = fieldPhase;

if nargin<3
    PulseOut.Gx = zeros(1,numSteps);
    PulseOut.Gy = zeros(1,numSteps);
    PulseOut.Gz = zeros(1,numSteps);
else
    switch (pulseSpatial)
        case 1,
            PulseOut.Gx = pulse1.Gx;
            PulseOut.Gy = pulse1.Gy;
            PulseOut.Gz = pulse1.Gz;
        case 2,
            PulseOut.Gx = pulse2.Gx;
            PulseOut.Gy = pulse2.Gy;
            PulseOut.Gz = pulse2.Gz;
        otherwise,
            error('pulseSpatial value %d should be 1 or 2', pulseSpatial);
    end
end