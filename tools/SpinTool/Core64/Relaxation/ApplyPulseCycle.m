function spins = ApplyPulseCycle(spins, pulse, pulsePhase, addCoeff, cycleType, B1)
% ApplyPulseCycle  Applies a shaped pulse with optional phase cycling
%   spins = ApplyPulseCycle(spins, pulse)  Applies the shaped pulse 
%   'pulse' to the input spin structure. 
%
%   spins = ApplyPulseCycle(spins, pulse, phaseCycle) Allows for phase
%   cycling. pulsePhase is a 1xN vector: N pulses are executed in parallel 
%   with phases pulsePhase(i), i=1,..,N, and the resulting magnetizations 
%   of the spins are added up. Phases are in degrees.
%
%   spins = ApplyPulseCycle(spins, pulse, phaseCycle, addCoeff)
%   Allows adding the different phase cycles with different weighting 
%   coefficients (including complex numbers!) in the 1xN vector addCoeff.
%
%   spins = ApplyPulseCycle(spins, pulse, phaseCycle, addCoeff, cycleType)
%   cycleType is a string which allows one to execute non-physical phase
%   cycles. Possible values (case insensitive):
%     'normal'      Phase cycle transformation R(phi) = Rz(phi)*R(0)*Rz(-phi)
%                   Here for example R(0)-R(pi/2)+R(pi)-R(3pi/2) yields the
%                   -1-->+1 and 1-->(-1) transitions ONLY.
%     'inverted'    Phase cycle transformation R(phi) = Rz(phi)*R(0)*Rz(+phi)
%                   Here for example R(0)+R(pi/2)+R(pi)+R(3pi/2) yields the
%                   -1-->+1, 1-->(-1) and 0-->0 transitions, so retains
%                   the longitudinal magnetization.
%
%   spins = ApplyPulseCycle(spins, pulse, phaseCycle, addCoeff, cycleType, B1)
%   B1 is a fixed value of B1. This can be used, for example, to apply hard
%   pulses which approximate adiabatic pulses by being B1 insensitive,
%   by simply fixing B1=1. Setting this to [] ignores the field.


RFphase = pulse.RFphase;

if nargin<3, pulsePhase = 0; end
numCycles = numel(pulsePhase);
if nargin<4, addCoeff = ones(1, numCycles)/numCycles; end
if nargin<5, cycleType = 'normal'; end
if nargin<6, B1 = []; end

numCycles = size(pulsePhase, 2);
numSpins = numel(spins);

% If user provided B1 value, we'll fix all spins to have that B1 value
if ~isempty(B1)
    tempB1 = zeros(1, numSpins);
    for idxSpin=1:numSpins
        tempB1(idxSpin) = spins(idxSpin).B1;
        spins(idxSpin).B1 = B1;
    end
end

if numCycles==1 % No phase cycle
    pulse.RFphase = pulse.RFphase + pulsePhase(1,1)/180*pi; 
    spins = ApplyPulseRelax(spins, pulse);
else % Phase cycling enabled
    spinsInitial = spins;
    for idxSpin=1:numSpins
        spins(idxSpin).M = [0;0;0];
    end
    switch lower(cycleType)
        case 'normal'  
            for idxCycle=1:numCycles
                curPhase = pulsePhase(idxCycle)/180*pi;
                curPulse = pulse;
                curPulse.RFphase = RFphase + curPhase;
                spinsOut = ApplyPulseRelax(spinsInitial, curPulse);
                for idxSpin=1:numSpins
                    spins(idxSpin).M = spins(idxSpin).M + addCoeff(idxCycle)*spinsOut(idxSpin).M;
                end
            end
        case 'inverted' 
            for idxCycle=1:numCycles
                curPhase = pulsePhase(idxCycle); % deg.
                spinsOut = ApplyZRotation(spinsInitial, curPhase);
                spinsOut = ApplyPulseRelax(spinsOut, pulse);
                spinsOut = ApplyZRotation(spinsOut, curPhase);
                for idxSpin=1:numSpins
                    spins(idxSpin).M = spins(idxSpin).M + addCoeff(idxCycle)*spinsOut(idxSpin).M;
                end
            end
        otherwise
            error('Unrecognized phase cycle %s', cycleType);
    end
end

% If user provided B1 value, we'll return all spins to their original B1 value
if ~isempty(B1)
    for idxSpin=1:numSpins
        spins(idxSpin).B1 = tempB1(idxSpin);
    end
end


