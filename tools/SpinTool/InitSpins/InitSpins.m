function spins = InitSpins(offsetVec, B1Vec, posVec, initMag, T1, T2, profileVec)
% INITSPINS  Initializes a spin structure
%
% Input:
%
% Name          Type         Units    Description      
% offsetVec     1xN double   kHz      Vector of chemical shifts
% B1Vec         1xM double   -        Scaling of B1
% posVec        1xP double   mm       Size of sample
% initMag       3x1          a.u.     Initial magnetization
% T1, T2        1x1 double   ms       Relaxation times
% profileVec    1xP double   -        Spatial profile of sample. Optional.
%                                     Modulates initMag spatially.

counter = 0;
if (nargin<7)
    profileVec = ones(1,numel(posVec));
end

% Initialize for speed
numTotalSpins = numel(offsetVec)*numel(B1Vec)*numel(posVec);
spins(numTotalSpins).r = [0; 0; 0];
spins(numTotalSpins).M  = initMag;
spins(numTotalSpins).cs = 0;  % in kHz!
spins(numTotalSpins).T1 = T1;
spins(numTotalSpins).T2 = T2;
spins(numTotalSpins).M0 = 1;
spins(numTotalSpins).B1 = 1;
spins(numTotalSpins).B0 = 0;
spins(numTotalSpins).RS = 1;


for idxOffset=1:numel(offsetVec)
    for idxB1=1:numel(B1Vec)
        for idxPos=1:numel(posVec)
            counter = counter + 1;
            spins(counter).r = [0; 0; posVec(idxPos)];
            spins(counter).M  = initMag.*profileVec;
            spins(counter).cs = offsetVec(idxOffset);  % in kHz!
            spins(counter).T1 = T1;
            spins(counter).T2 = T2;
            spins(counter).M0 = 1;
            spins(counter).B1 = B1Vec(idxB1);
            spins(counter).B0 = 0;
            spins(counter).RS = 1;
        end
    end;
end;
