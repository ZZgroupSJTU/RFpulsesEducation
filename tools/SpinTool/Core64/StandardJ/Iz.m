function M = Iz(spin)
% Returns the normalized Iz spin operator.
%   M = Iz  returns the 2x2 z-spin operator for a spin-1/2. The operator
%   is normalized according to the Pauli matrices' commutation relations:
%     [Ix, Iy] = 2iIz
%     [Iy, Iz] = 2iIx
%     [Iz, Ix] = 2iIy
%   M = Iz(spin)  Returns the (normalized) Pauli matrix for a given spin.
%   Allowed values: 1/2, 1, 3/2, 2.

if nargin<1, spin = 1/2; end

switch spin
    case 1/2
        M = [1 0; 0 -1]/2;
    case 1
        M = [1 0  0;
             0 0  0;
             0 0 -1];
    case 3/2
        M = [3/2   0    0    0;
               0 1/2    0    0;
               0   0 -1/2    0;
               0   0    0 -3/2];
    case 2
        M = [2  0  0  0  0;
             0  1  0  0  0;
             0  0  0  0  0;
             0  0  0 -1  0;
             0  0  0  0 -2];
    otherwise
        error('Can only generate spin matrices up to spin 2.');
end


