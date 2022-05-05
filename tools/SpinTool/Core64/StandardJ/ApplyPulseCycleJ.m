function spins = ApplyPulseCycleJ(spins, pulse, affectedNuclei, pulsePhase, addCoeff, freqRange)
% Applies a shaped pulse with optional phase cycling.
%   spins = ApplyPulseCycle(spins, pulse)  Applies the shaped pulse 
%   'pulse' to the input spin structure. 
%
%   spins = ApplyPulseCycle(spins, pulse, pulsePhase) Allows for phase
%   cycling. pulsePhase is a 1xN vector: N pulses are executed in parallel 
%   with phases pulsePhase(i), i=1,..,N, and the resulting magnetizations 
%   of the spins are added up. Phases are in degrees.
%
%   spins = ApplyPulseCycle(spins, pulse, pulsePhase, addCoeff)
%   Allows adding the different phase cycles with different weighting 
%   coefficients (including complex numbers!) in the 1xN vector addCoeff.
%
%   spins = ApplyPulseCycle(spins, pulse, pulsePhase, addCoeff,freqRange)
%   freqRange is an optional 1x2 vector of min. and max. ppm affected by
%   the pulse. This is a crude but convenient way of simulating selective
%   pulses.


RFphase = pulse.RFphase;
isAcquire = 0;

if nargin<3, affectedNuclei = []; end
if nargin<4, pulsePhase = []; end
if nargin<5, addCoeff = []; end
if nargin<6, freqRange = []; end

if isempty(pulsePhase), pulsePhase = 0; end
numCycles = numel(pulsePhase);
if isempty(addCoeff), addCoeff = ones(1, numCycles)/numCycles; end
numCycles = size(pulsePhase, 2);
numMolecules = numel(spins.molecule);

if numCycles==1 % No phase cycle
    pulse.RFphase = pulse.RFphase + pulsePhase(1,1)/180*pi; 
    [~, spins] = PropagateJ(spins, pulse, isAcquire, affectedNuclei, freqRange);
else % Phase cycle
    spinsInitial = spins;
    % Let's zero out rho in preparation for adding up the results of the different cycles.
    for idxMolecule=1:numMolecules
        zeroRho = [];
        numSpins = numel(spins.molecule(idxMolecule).spin);
        if iscell(spins.molecule(idxMolecule).csVec)
            numSubSystems = numel(spins.molecule(idxMolecule).csVec);
            for idxSys=1:numSubSystems
                numNuclei = numel(spins.molecule(idxMolecule).csVec{idxSys});
                zeroRho{idxSys} = zeros(2^numNuclei);
            end
        else
            numNuclei = numel(spins.molecule(idxMolecule).csVec);
            zeroRho = zeros(2^numNuclei);
        end
        for idxSpin=1:numSpins
            spins.molecule(idxMolecule).spin(idxSpin).rho = zeroRho;
        end
    end
    % Perform each cycle and add the result to the density matrix of each spin (of each molecule).
    for idxCycle=1:numCycles
        curPhase = pulsePhase(idxCycle)/180*pi;
        curPulse = pulse;
        curPulse.RFphase = RFphase + curPhase;
        [~, spinsOut] = PropagateJ(spinsInitial, curPulse, isAcquire, affectedNuclei, freqRange); 
        for idxMolecule=1:numMolecules
            numSpins = numel(spins.molecule(idxMolecule).spin);
            if iscell(spins.molecule(idxMolecule).csVec)
                numSubSystems = numel(spins.molecule(idxMolecule).csVec);
                for idxSys=1:numSubSystems
                    for idxSpin=1:numSpins
                        spins.molecule(idxMolecule).spin(idxSpin).rho{idxSys} = ...
                            spins.molecule(idxMolecule).spin(idxSpin).rho{idxSys} + ...
                            addCoeff(idxCycle)*spinsOut.molecule(idxMolecule).spin(idxSpin).rho{idxSys};
                    end
                end
            else
                for idxSpin=1:numSpins
                    spins.molecule(idxMolecule).spin(idxSpin).rho = spins.molecule(idxMolecule).spin(idxSpin).rho + addCoeff(idxCycle)*spinsOut.molecule(idxMolecule).spin(idxSpin).rho;
                end
            end
        end
    end
end



