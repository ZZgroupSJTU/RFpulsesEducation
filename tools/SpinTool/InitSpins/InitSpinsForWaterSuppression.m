function spins = InitSpinsForWaterSuppression(T1Vec, B1Vec, T2, chemShift)
%   spins = InitSpinsForWaterSuppression(T1Vec, B1Vec, T2, chemShift)
%
%   The current function creates a spin structure with no J coupling, with
%   a specified range of T1 and B1-inhomogeneity values 
%
%   Inputs 
%
%   Name              Type         Units   Description      
%   T1Vec             1xN double   ms      Vector of T1 values
%   B1Vec             1xM double   -       A vector of B1 scaling values
%   T2                double       ms      Transverse relaxation
%   chemShift         double       kHz     chemical shift offset of Water

% Populate spin structure
counter = 0;
for idxT1=1:numel(T1Vec)
    for idxB1=1:numel(B1Vec)
        counter = counter + 1;
        spins(counter).r  = [0; 0; 0];
        spins(counter).M  = [0; 0; 1];
        spins(counter).cs = chemShift;  % in kHz!
        spins(counter).T1 = T1Vec(idxT1); % ms
        spins(counter).T2 = T2; % ms
        spins(counter).M0 = 1; % a.u.
        spins(counter).B1 = B1Vec(idxB1); % Scales RF
        spins(counter).B0 = 0; % Offset, in kHz
        spins(counter).RS  = 1; % Receiver sensitivity
    end;
end;
