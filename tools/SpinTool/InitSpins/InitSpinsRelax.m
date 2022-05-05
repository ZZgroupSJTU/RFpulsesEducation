function spins = InitSpinsRelax(chemShiftVec, ...
                                numSpinsPerShift, ...
                                sampleSize, ...
                                initialMag, ...
                                T1, ...
                                T2, ...
                                eqMag, ...
                                sampleOffset, ...
                                B1, ...
                                spatialAxis)
%   spins = InitSpinsRelax(chemShiftVec,numSpinsPerShift,sampleSize,initialMag,T1,T2,eqMag, sampleOffset, B1, spatialAxis)
%
%   The current function creates a 1D spin structure with no J coupling, to
%   which pulses, delays, etc. can be acurSpinlied using AcurSpinlyPulseRelax, etc.
%
%   Inputs 
%
%   Name              Type         Units   Description      
%   chemShiftVec      1xN double   kHz     Vector of chemical shifts
%   numSpinsPerShift  double       -       # of spins per chemical shift
%   sampleSize        double       mm      Size of sample
%   initialMag        3x1 double   -       Initial magnetization
%   T1                double       ms      Longitudinal relaxation
%   T2                double       ms      Transverse relaxation
%   eqMag             double       -       Equilibrium magnetization (M0)
%   sampleOffset      double       mm      Offset of the sample's center.
%                                          Optional (set to 0 if omitted)
%   B1                double       -       B1 scaling factor
%   spatialAxis       char         -       'x', 'y', or 'z'. Specifies 
%                                          along which of the axes the
%                                          spins will be created
%
%   Example:
%   spins = InitSpinsRelax(0, 100, 100, [0; 0; 1], 1000, 200, 1);

numShifts = length(chemShiftVec);
dz = sampleSize/numSpinsPerShift;
if (numSpinsPerShift>1)
    % posVec = [-sampleSize/2:dz:sampleSize/2-dz];
    posVec = linspace(-sampleSize/2, sampleSize/2, numSpinsPerShift);
else
    posVec = 0;
end

if nargin<10
    spatialAxis = 'z';
else
    if isempty(spatialAxis)
        spatialAxis = 'z';
    end
end
spatialAxis = lower(spatialAxis);

if (nargin<9)
    B1 = 1;
end

if (nargin<8)
    sampleOffset = 0;
end

if (nargin<7)
    eqMag = 1;
end

if (nargin<6)
    T2 = 1e6;
end

if (nargin<5)
    T1 = 1e6;
end

% Pre-allocate for speed
numTotalSpins = numel(chemShiftVec)*numSpinsPerShift;
spins(numTotalSpins).r = [0; 0; 0];
spins(numTotalSpins).M  = initialMag;
spins(numTotalSpins).cs = chemShiftVec(1);  % in kHz!
spins(numTotalSpins).T1 = T1;
spins(numTotalSpins).T2 = T2;
spins(numTotalSpins).M0 = eqMag;
spins(numTotalSpins).B1 = B1;
spins(numTotalSpins).B0 = 0;
spins(numTotalSpins).RS = 1;

% Populate spin structure
counter = 0;
for curShift=1:numShifts
    for curSpin=1:numSpinsPerShift
        counter = counter + 1;
        switch spatialAxis
            case 'x'
                spins(counter).r  = [posVec(curSpin) + sampleOffset; 0; 0];
            case 'y'
                spins(counter).r  = [0; posVec(curSpin) + sampleOffset; 0];
            case 'z'
                spins(counter).r  = [0; 0; posVec(curSpin) + sampleOffset];
            otherwise
                error('Unknown spatialAxis.');
        end
        spins(counter).M  = initialMag;
        spins(counter).cs = chemShiftVec(curShift);  % in kHz!
        spins(counter).T1 = T1; % ms
        spins(counter).T2 = T2; % ms
        spins(counter).M0 = eqMag; % a.u.
        spins(counter).B1 = B1; % Scales RF
        spins(counter).B0 = 0; % Offset, in kHz
        spins(counter).RS  = 1; % Receiver sensitivity
    end;
end;
