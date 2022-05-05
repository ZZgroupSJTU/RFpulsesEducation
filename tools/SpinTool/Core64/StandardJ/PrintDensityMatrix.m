function PrintDensityMatrix(DM, epsilon)
% PrintDensityMatrix(DM, epsilon)
%
%    Prints the density matrix of a 2-spin system using product operator
%    formalism. The optional epsilon parameter (which defaults to 0.001
%    if omitted) chooses a weighting cutoff threshold, below which 
%    product operator elements will not be displayed.

% Note to self:
% Don't forget Ix, Iy, Iz, Sx, Sy, Sz used in the literature have a
% factor of 1/2 included, as opposed to the Pauli Matrices used in our
% simulations.

s='';
if (nargin<2)
    epsilon = 0.001;
end

x = trace(DM*kron(Ix,I));
if (abs(x)>epsilon)
    s = sprintf('%s%.3f Ix + ',s,x);
end

x = trace(DM*kron(Iy,I));
if (abs(x)>epsilon)
    s = sprintf('%s%.3f Iy + ',s,x);
end

x = trace(DM*kron(Iz,I));
if (abs(x)>epsilon)
    s = sprintf('%s%.3f Iz + ',s,x);
end

% ------------------------------------------------------------------------

x = trace(DM*kron(I,Ix));
if (abs(x)>epsilon)
    s = sprintf('%s%.3f Sx + ',s,x);
end

x = trace(DM*kron(I,Iy));
if (abs(x)>epsilon)
    s = sprintf('%s%.3f Sy + ',s,x);
end

x = trace(DM*kron(I,Iz));
if (abs(x)>epsilon)
    s = sprintf('%s%.3f Sz + ',s,x);
end

% ------------------------------------------------------------------------

x = trace(DM*kron(Ix,Ix))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IxSx + ',s,x);
end

x = trace(DM*kron(Ix,Iy))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IxSy + ',s,x);
end

x = trace(DM*kron(Ix,Iz))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IxSz + ',s,x);
end

% ------------------------------------------------------------------------

x = trace(DM*kron(Iy,Ix))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IySx + ',s,x);
end

x = trace(DM*kron(Iy,Iy))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IySy + ',s,x);
end

x = trace(DM*kron(Iy,Iz))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IySz + ',s,x);
end

% ------------------------------------------------------------------------

x = trace(DM*kron(Iz,Ix))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IzSx + ',s,x);
end

x = trace(DM*kron(Iz,Iy))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IzSy + ',s,x);
end

x = trace(DM*kron(Iz,Iz))*4;
if (abs(x)>epsilon)
    s = sprintf('%s%.3f IzSz + ',s,x);
end

% ------------------------------------------------------------------------

if isempty(s)
    % One can get zero if the density matrix is non hermitian (and
    % there invalid)
    fprintf('Zero / invalid density matrix\n',s);
else
    s = s(1:end-3); % Remove the ' + ' from the end of the string
    fprintf('%s\n',s);
end