function M=SumSpins(spins, isNormalize)
% M=SumSpins(spins, isNormalize)
%
% Sums the magnetization vectors over all spins in the ensemble.
% If isNormalize is set to true, the result will be normalized by the
% number of spins in the ensemble. 

if nargin<2, isNormalize = false; end

M = [0; 0; 0];
for idx=1:numel(spins)
    M = M + spins(idx).M;
end
if isNormalize, M = M/numel(spins); end