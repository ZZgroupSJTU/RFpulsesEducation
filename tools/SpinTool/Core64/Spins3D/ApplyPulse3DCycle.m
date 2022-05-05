function spins = ApplyPulse3DCycle(spins, pulse, pulsePhase, addCoeff, cycleType, B1)
% ApplyPulse3DCycle  Applies a shaped pulse with optional phase cycling.
%   Consult ApplyPulseCycle for possible syntax and optional inputs to
%   this function. This function does exactly the same, only operating 
%   on 3D spin structures (as created by InitSpins3D).

RFphase = pulse.RFphase;

if nargin<3, pulsePhase = 0; end
numCycles = numel(pulsePhase);
if nargin<4, addCoeff = ones(1, numCycles)/numCycles; end
if nargin<5, cycleType = 'normal'; end
if nargin<6, B1 = []; end

% If user provided B1 value, we'll fix all spins to have that B1 value
if ~isempty(B1)
    tempB1 = spins.B1;
    spins.B1 = B1;
end


numCycles = size(pulsePhase, 2);
numSpins = numel(spins);
if numCycles==1 % No phase cycle
    pulse.RFphase = pulse.RFphase + pulsePhase(1,1)/180*pi; 
    spins = ApplyPulse3D(spins, pulse);
else % Phase cycling enabled
    spinsInitial = spins;
    for idxSpin=1:numSpins
        spins.M = zeros(size(spins.M));
    end
    switch lower(cycleType)
        case 'normal'  
            for idxCycle=1:numCycles
                curPhase = pulsePhase(idxCycle)/180*pi;
                curPulse = pulse;
                curPulse.RFphase = RFphase + curPhase;
                spinsOut = ApplyPulse3D(spinsInitial, curPulse);
                spins.M = spins.M + addCoeff(idxCycle)*spinsOut.M;
            end
        case 'inverted' 
            for idxCycle=1:numCycles
                curPhase = pulsePhase(idxCycle); % deg.
                spinsOut = ApplyZRotation3D(spinsInitial, curPhase);
                spinsOut = ApplyPulse3D(spinsOut, pulse);
                spinsOut = ApplyZRotation3D(spinsOut, curPhase);
                spins.M = spins.M + addCoeff(idxCycle)*spinsOut.M;
            end
        otherwise
            error('Unrecognized phase cycle %s', cycleType);
    end
end

% If user provided B1 value, we'll return the B1 of the spins to its original value
if ~isempty(B1)
    spins.B1 = tempB1;
end


