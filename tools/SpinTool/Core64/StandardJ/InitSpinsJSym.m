function spins = InitSpinsJSym(csVec, csCenter, JMatrix, B0, relAmplitudes, sampleSize, numMolecules)
% SYNTAX: spins = InitSpinsJSym(csVec, JMatrix, relAmplitudes, sampleSize, numMolecules)
%
% Input Variables
% Variable Name  Size         Units   Description
% rho            (2^M)x(2^M)  -       Density matrix of system
% csVec          1xM          ppm     Chemical shift of each Hydrogen atom
% csCenter       1x1          ppm     Center frequency of Tx/Rx
% JMatrix        MxM          Hz      A (symmetric) matrix of J-couplings
% B0             1x1          Tesla   Main field
% relAmplitudes  1xM          -       Relative amplitudes of each nucleus
% sampleSize     1x1          mm      Size of sample
% numMolecules   1x1          -       Number of molecules, dispersed 
%                                     evenly within the sample along z.
%
% Example:
%
% To initialize a GABA molecule (A2M2X2 system), use:
% csVec = [3.01 3.01 2.28 2.28 1.89 1.89]
% csCenter = 4.7; % Makes water in the center
% J = 7; % Hz (rough estimate for GABA)
% JMatrix = [0 0 0 0 J J;
%            0 0 0 0 J J;
%            0 0 0 0 J J;
%            0 0 0 0 J J;
%            J J J J 0 0;
%            J J J J 0 0];
% relAmplitudes = [1 1 1 1 1 1]; % Equal amplitudes for all nuclei
% B0 = 3; % Tesla
% spins = InitSpinsJ(csVec, csCenter, JMatrix, B0, relAmplitudes);

zAxis = linspace(-sampleSize/2, sampleSize/2, numMolecules); % mm
numSpins = numel(csVec);

% Create equilibrium density matrix
rho0 = zeros(2^numSpins, 2^numSpins);
for idx=1:numSpins
    rho0 = rho0 + relAmplitudes(idx)*IzN(idx, numSpins);
end

spins(numMolecules).r = [0; 0; 0];
for idxMolecule=1:numMolecules
    spins(idxMolecule).csVec = csVec;
    spins(idxMolecule).csCenter = csCenter;
    spins(idxMolecule).JMatrix = JMatrix;
    spins(idxMolecule).r = [0; 0; zAxis(idxMolecule)];
    spins(idxMolecule).B0 = B0; 
    spins(idxMolecule).rho = rho0;
end

