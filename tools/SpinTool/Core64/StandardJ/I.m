function M = I(spin)
% Identity matrix of given spin.
%   M = I  2x2 identity matrix for a spin-1/2;
%
%   M = I(spin)  Identity matrix for a given spin. For example, I(3/2)
%   is the identity matrix for a spin-3/2, which is a 4x4 identity matrix.

if nargin<1, spin = 1/2; end
numDim = round(spin*2+1);
M = eye(numDim); 