function [spins, fid] = ApplyPulseExchangeCycle(spins, pulse, phaseCycle, addCoeff, pulseOffset)
% ApplyPulseExchangeCycle  Applies a shaped pulse with optional phase cycling
%   spins = ApplyPulseExchangeCycle(spins, pulse) Applies the shaped pulse 
%   'pulse' to the input spin structure. This basically has no cycling
%   and just applies the pulse once "as-is". If you're doing this, you're
%   probably better off using ApplyPulseExchange directly.
%
%   spins = ApplyPulseExchangeCycle(spins, pulse, phaseCycle) Allows for phase
%   cycling. phaseCycle is a 1xN vector: N pulses are executed in parallel 
%   with phases phaseCycle(i), i=1,..,N, and the resulting magnetizations 
%   of the spins are normalized by N and added up. Phases are in degrees.
%
%   spins = ApplyPulseExchangeCycle(spins, pulse, phaseCycle, addCoeff)
%   Allows adding the different phase cycles with different weighting 
%   coefficients (including complex numbers!) in the 1xN vector addCoeff.
% 
%   spins = ApplyPulseExchangeCycle(spins, pulse, phaseCycle, addCoeff, pulseOffset)
%   Allows adding a constant offset (in kHz) to the pulse. If omitted, offset
%   defaults to 0.


if nargin<3, phaseCycle = 0; end
numCycles = numel(phaseCycle);
if nargin<4, addCoeff = ones(1, numCycles)/numCycles; end
if nargin<5, pulseOffset = 0; end

numCycles = numel(phaseCycle);
numSpins = numel(spins);

if numCycles==1 
    % No phase cycle (single pulse)
    [spins, fid] = ApplyPulseExchange(spins, pulse, pulseOffset, phaseCycle); 
else
    % Phase cycle through pulses.
    % Initialize variables.
    numSteps = numel(pulse.RFamp);
    fid = zeros(1, numSteps);
    spinsInitial = spins;
    for idxSpin=1:numSpins
        spins(idxSpin).M = 0.*spins(idxSpin).M;
    end
    % Apply the cycle
    for idxCycle=1:numCycles
        [spinsOut, fidOut] = ApplyPulseExchange(spinsInitial, pulse, pulseOffset, phaseCycle(idxCycle));
        for idxSpin=1:numSpins
            spins(idxSpin).M = spins(idxSpin).M + addCoeff(idxCycle)*spinsOut(idxSpin).M;
        end
        fid = fid + addCoeff(idxCycle)*fidOut;
    end
end



