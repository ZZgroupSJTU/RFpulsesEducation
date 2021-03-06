function spins = InitSpinsForWaterSuppressionB0B1(B0Vec, B1Vec, T1, T2)
%   spins = InitSpinsForWaterSuppression(T1Vec, B1Vec, T2)
%
%   The current function creates a spin structure with no J coupling, with
%   a specified range of B0 and B1-inhomogeneity values 
%
%   Inputs 
%   Name              Type         Units   Description      
%   B0Vec             1xN double   kHz     Vector of B0 ("chem. shift") values
%   B1Vec             1xM double   -       A vector of B1 scaling values
%   T1                double       ms      Longitudinal relaxation
%   T2                double       ms      Transverse relaxation

% Populate spin structure
counter = 0;
for idxB0=1:numel(B0Vec)
    for idxB1=1:numel(B1Vec)
        counter = counter + 1;
        spins(counter).r  = [0; 0; 0];
        spins(counter).M  = [0; 0; 1];
        spins(counter).cs = B0Vec(idxB0);  % in kHz!
        spins(counter).T1 = T1; % ms
        spins(counter).T2 = T2; % ms
        spins(counter).M0 = 1; % a.u.
        spins(counter).B1 = B1Vec(idxB1); % Scales RF
        spins(counter).B0 = 0; % Offset, in kHz
        spins(counter).RS  = 1; % Receiver sensitivity
    end;
end;
