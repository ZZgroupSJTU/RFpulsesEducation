function spins = ApplyInhoToSpins(spins, coeff, inhoAxis, applyToWhat)
% ApplyInhoToSpins  Applies polynomial B0 or B1 inhomogeneity to spins.
%   Applies a polynomial B1 inhomogeneity to the spins, as a function
%   of position. The current function only applies the inhomogeneity along
%   the specified axis ('x', 'y', 'z'). Note the coefficients in coeff
%   must be specified from HIGHEST order to LOWEST (zeroth).
%
% Input Variables
% Variable Name      Description
% spins              Input spin structure
% coeff              Vector containing polynomial coefficients of
%                    inhomogeneity (in units of Hz/mm). The coefficients 
%                    must be specified from HIGHEST order to LOWEST (0th)
% inhoAxis           'x', 'y', 'z' (case insensitive)
% applyToWhat        'b0', 'b1' (case insensitive). Default: B0.

numSpins = numel(spins);
if nargin<4, applyToWhat = 'B0'; end

% Extract spins' positions
for idx=1:numSpins
    r(idx, :) = spins(idx).r;
end

% Calculate inhomogeneity for each spin
switch lower(inhoAxis)
    case 'x'
        N = numel(coeff); 
        V = bsxfun(@power, r(1,:).', 0:N);
        inho = V*coeff;
    case 'y'
        N = numel(coeff); 
        V = bsxfun(@power, r(2,:).', 0:N);
        inho = V*coeff;
    case 'z'
        N = numel(coeff); 
        V = bsxfun(@power, r(3,:).', 0:N);
        inho = V*coeff;
    otherwise
        error('Unrecognized inhoAxis = %s', inhoAxis);
end

% Assign the inhomogeneity to the spins
switch lower(applyToWhat)
    case 'b0'
        for idx=1:numSpins
            spins(idx).B0 = inho(idx);
        end
    case 'b1'
        for idx=1:numSpins
            spins(idx).B1 = inho(idx);
        end
    otherwise
        error('Error in ApplyB1InhoToSpins: applyToWhat not B0,B1! Aborting.');
end
