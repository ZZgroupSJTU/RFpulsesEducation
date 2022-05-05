function spins = SpinsJAddSpin(spins, molecule, B0, relAmplitude, isSecular, linewidth)
% SYNTAX: spins = SpinsJAddSpin(spins, molecule, isSecular, linewidth)
%
% Adds a specific molecule type to the spin structure.
%
% Input Variables
% spins        The input spin structure
% molecule     The input molecule. A string. One of the following (case
%              insensitive):
%              'Spin Half': On resonance spin-1/2.
%              'Simple J': Two spin-1/2s, at +-1 ppm with J=10 Hz.
%              'NAA'
%              'GABA'
%              'Lac'
%              'mI'
%              'Cho'
%              'Cr'
%              'Cho singlet'
%              'Cr singlet'
%              'NAA singlet'
%              'Glu'
%              'Gln'
%
% Output Variables
% spins        The input spin structure with the appended molecule.

numSpins = numel(spins);
if nargin<3, B0 = 3; end
if nargin<4, relAmplitude = 1; end
if nargin<5, isSecular = 0; end
if nargin<6, linewidth = 5; end

switch lower(molecule)
    case 'spin half'
        csVec = 0;
        csCenter = 0;
        JMatrix = 0;
        relAmplitudes = relAmplitude;
        sampleSize = 0;
        numMolecules = 1;
        nuclei = {'1H'};
        B1 = 1;
        spins(numSpins+1) = InitSpinsJ(csVec, csCenter, JMatrix, B0, relAmplitudes, sampleSize, numMolecules, nuclei, B1, isSecular, linewidth);
    case 'simple j'
        csVec = [-0.5 0.5]; % ppm
        csCenter = 0; % ppm
        J = 10; % Hz
        JMatrix = [0 J;
                   J 0]; % Hz
        relAmplitudes = [1 1].*relAmplitude;
        sampleSize = 0;
        numMolecules = 1;
        nuclei = {'1H', '1H'};
        B1 = 1; % Scaling
        spins(numSpins+1) = InitSpinsJ(csVec, csCenter, JMatrix, B0, relAmplitudes, sampleSize, numMolecules, nuclei, B1, isSecular, linewidth);
    case 'naa singlet'
        csVec = 0;
        csCenter = 0;
        JMatrix = 0;
        relAmplitudes = relAmplitude;
        sampleSize = 0;
        numMolecules = 1;
        nuclei = {'1H'};
        B1 = 1;
        spins(numSpins+1) = InitSpinsJ(csVec, csCenter, JMatrix, B0, relAmplitudes, sampleSize, numMolecules, nuclei, B1, isSecular, linewidth);
    otherwise
        error('Unidentified molecule type in SpinsJAddSpin.');
end