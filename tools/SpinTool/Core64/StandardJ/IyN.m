function M = IyN(idx, numNuclei, spin)
% Creates an Iy operator for a given spin.
%   M = IyN(idx, numNuclei) Create the Iy spin operator for the idx spin
%   in a spin system having numNuclei spin-half spins. The resulting
%   operator is a 2^numNuclei * 2^numNuclei matrix.
%
%   M = IyN(idx, numNuclei, spin) The optional vector input spin describes
%   the spin of each of the nuclei in the system. For example, if you have
%   a system with 1H, 14N and 1H nuclei, then spin=[1/2, 1, 1/2]. 

if nargin<3
    spin = ones(1,numNuclei)*1/2;
end

if (idx==1)
    M = Iy(spin(1));
else
    M = I(spin(1));
end

for k=2:numNuclei
    if (k==idx)
        M = kron(M,Iy(spin(k)));
    else
        M = kron(M,I(spin(k)));
    end
end

