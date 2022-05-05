function CSA = CalculateCSA(BW, cs)
% CSA = CalculateCSA(BW, cs)
%
% Calculates the chemical shift artifact associated with a pulse of a 
% certain bandwidth and a chemical shift cs. Both cs and BW must have the
% same units.
%
% The chemical shift artifact is the shift dz, in mm, given by
%
% [1] gamma*G*dz = v0
%
% where v0 is the chemical shift and G is the gradient needed to 
% excite a sample of length L, gamma*G*L=BW. Substituting this in the
% above, we obtain
%
% [2] CSA = dz/L = v0/BW
% 
% where CSA is the fractional shift, dz/L, independent of the size of 
% the length L excited. To calculate the shift for a slice of width L,
% just calculate CalculateCSA(BW,cs)*L.

CSA = cs/BW;
