function spinsOut = PurgePerfect(spinsIn, spoilingDuration)
% Description: artificially sets to 0 the transverse components of the
% mangetization.

spinsOut=spinsIn;
numSpins = length(spinsIn);
for idx=1:numSpins
    spinsOut(idx).M(1) = 0;
    spinsOut(idx).M(2) = 0;
end

if (nargin>1)
    spinsOut = DelayRelax(spinsOut, spoilingDuration);
end