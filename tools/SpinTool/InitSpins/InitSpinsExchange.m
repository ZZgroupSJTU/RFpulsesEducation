function spins = InitSpinsExchange(chemShiftVec, ...
                                   exchangeMatrix, ...
                                   initialMag, ...
                                   T1, ...
                                   T2, ...
                                   eqMag, ...
                                   B1Scaling, ...
                                   receiverSensitivity)
% InitSpinsExchange
%   The current function creates a 1D spin structure with no J coupling, and 
%   with chemical exchange.
%
%   Inputs 
%
%   Name              Type         Units   Description      
%   chemShiftVec      1xN double   kHz     Vector of chemical shifts
%   exchangeMatrix    NxN matrix   kHz     Exchange matrix between sites
%   initialMag        3Nx1 double  -       Initial magnetization
%   T1                1xN double   ms      Longitudinal relaxation
%   T2                1XN double   ms      Transverse relaxation
%   eqMag             1xN double   -       Equilibrium magnetization (M0)
%   B1                1xN double   -       Transmitter sensitivity (B1+)
%   RS                1xN double   -       Receiver sensitivity (B1-)
%
%   where N = number of exchange sites (chemical shifts).
%   T1, T2, B1, RS, eqMag, can also be just numbers, in which case they
%   are assumed to be the same for all exchanging sites.
%   All spins are assumed to be at x=y=z=0. The position of each spin
%   can be accessed via the r field, e.g.: spins(j).r = [3;2;5]
%
%   Example:
%   k = 1;  % Symmetric exchange at rate of 1 kHz
%   exchangeMatrix = [-k  k; 
%                      k -k];
%   chemShiftVec = [-0.5 0.5];  % kHz
%   T1 = 1000; % ms
%   T2 = 200; % ms
%   eqMag = [1; 2];           % Normalized. Both spins have same M0
%   B1Scaling = 1;            % No B1+ inhomogeneity
%   receiverSensitivity = 1;  % No B1- inhomogeneity
%   initialMag = [0; 0; 1; 0; 0; 1];  % Both spins start with equal magnitude from z-axis
%   spins = InitSpinsExchange(chemShiftVec, exchangeMatrix, initialMag, T1, T2, eqMag, B1Scaling, receiverSensitivity);

if nargin<3, initialMag = [0; 0; 1]; end
if nargin<4, T1 = 1e9; end
if nargin<5, T2 = 1e9; end
if nargin<6, eqMag = 1; end
if nargin<7, B1Scaling = 1; end
if nargin<8, receiverSensitivity = 1; end

numSites = size(exchangeMatrix, 1);

spins(1).r = [0; 0; 0];
if numel(initialMag)==3
    spins(1).M = repmat(initialMag, numSites, 1);
else
    spins(1).M = initialMag;
end
spins(1).cs = chemShiftVec;
spins(1).T1 = T1;
spins(1).T2 = T2;
spins(1).M0 = eqMag;
spins(1).B1 = B1Scaling;
spins(1).RS = receiverSensitivity;
spins(1).K = exchangeMatrix;

