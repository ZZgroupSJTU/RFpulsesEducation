function M = Ix(spin)
% Returns the normalized Ix spin operator.
%   M = Ix  returns the 2x2 x-spin operator for a spin-1/2. The operator
%   is normalized according to the Pauli matrices' commutation relations:
%     [Ix, Iy] = 2iIz
%     [Iy, Iz] = 2iIx
%     [Iz, Ix] = 2iIy
%   M = Ix(spin)  Returns the (normalized) Pauli matrix for a given spin.
%   Allowed values: 1/2, 1, 3/2, 2.

if nargin<1, spin = 1/2; end

switch spin
    case 1/2
        M = [0 1; 1 0]/2;
    case 1
        M = [0 1 0;
             1 0 1;
             0 1 0]*1/sqrt(2);
    case 3/2
        a = sqrt(3);
        M = [ 0  a  0 0;
              a  0  2 0;
              0  2  0 a;
              0  0  a 0]*1/2;
    case 2
        a = sqrt(6);
        M = 1/2*[ 0  2  0  0  0;
                  2  0  a  0  0;
                  0  a  0  a  0;
                  0  0  a  0  2;
                  0  0  0  2  0];
    otherwise
        error('Can only generate spin matrices up to spin 2.');
end





