function spins = AddSpinsExchange(spins, ...
                                  chemShiftVec, ...
                                  exchangeMatrix, ...
                                  initialMag, ...
                                  T1, ...
                                  T2, ...
                                  eqMag, ...
                                  B1Scaling, ...
                                  receiverSensitivity)
% spins = InitSpinsExchange(spins, ...
%                           chemShiftVec, ...
%                           exchangeMatrix, ...
%                           initialMag, ...
%                           T1, ...
%                           T2, ...
%                           eqMag, ...
%                           B1Scaling, ...
%                           receiverSensitivity)
%
%    The current function creates a 1D spin structure with no J coupling, to
%    which pulses, delays, etc. can be applied. Only a single "spin system"
%    (with possibly exchanging chemical sites) is created.
%
%    Inputs 
%
%    Name              Type         Units   Description      
%    spins             -            -       Input spin structure (use [] if none exists yet)
%    chemShiftVec      1xN double   kHz     Vector of chemical shifts
%    exchangeMatrix    NxN double   kHz     Exchange reaction rates
%    sampleSize        double       mm      Size of sample
%    initialMag        3x1 double   -       Initial magnetization
%    T1                double       ms      Longitudinal relaxation
%    T2                double       ms      Transverse relaxation
%    eqMag             double       -       Equilibrium magnetization (M0)
%    B1                double       -       B1 scaling factor

if nargin<4, initialMag = [0; 0; 1]; end
if nargin<5, T1 = 1e9; end
if nargin<6, T2 = 1e9; end
if nargin<7, eqMag = 1; end
if nargin<8, B1Scaling = 1; end
if nargin<9, receiverSensitivity = 1; end

numSites = size(exchangeMatrix, 1);

spinsSystem.r = [0; 0; 0];
if numel(initialMag)==3
    spinsSystem.M = repmat(initialMag, numSites, 1);
else
    spinsSystem.M = initialMag;
end
spinsSystem.cs = chemShiftVec;
spinsSystem.T1 = T1;
spinsSystem.T2 = T2;
spinsSystem.M0 = eqMag;
spinsSystem.B1 = B1Scaling;
spinsSystem.RS = receiverSensitivity;
spinsSystem.K = exchangeMatrix;

spins = [spins, spinsSystem];