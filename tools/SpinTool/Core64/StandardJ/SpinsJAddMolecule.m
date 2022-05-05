function spins = SpinsJAddMolecule(spins, molecule, concentration, r)
% Adds a specific molecule type to the spin structure. 
%   spins = SpinsJAddMolecule(spins, molecule) spins is either an existing
%   J-coupled spin structure, or an empty set [] (in which case a new 
%   default spin structure will be created). molecule is a string indicating
%   the molecule type, e.g. 'cho' or 'NAA' (case insensitive). For a complete
%   list of available molecules, type in SpinsJAddMolecule with no further
%   input. 
%
%   spins = SpinsJAddMolecule(spins, molecule, concentration) Allows the
%   user to specify a concentration. If not given or set to [], the 
%   molecule's concentration will be assumed based on standard in-vivo
%   values in humans.
%
%   spins = SpinsJAddMolecule(spins, molecule, concentration, r) Allows
%   the user to specify the position of the newly added molecule (mm). r is
%   a 3xN vector, where N = number of spins, and r(:,j) is the position of
%   the 1<=j<=N spin. If omitted or set to [], then r=[0;0;0] by default
%   (i.e. single spin at the origin).
%
%   For references regarding chemical shifts and J-coupling constants of
%   popular metabolites, consult:
%     [1] Govind et al, NMR Biomed 28:923-924 (2015)
%     [2] Govindaraju et al, NMR Biomed 13:129-153 (2000)

[~, moleculeName, moleculeAbbrev, humanConcentration] = SpinsJGetMolecule(molecule);

% If default values are not provided 
if nargin<3, concentration = humanConcentration; end
if nargin<4, r = [0;0;0]; end

% The relative amplitude of each molecule here is used to correct for 
% cases in which several protons are grouped together. For example, if 
% we create an NAA singlet which is derived from a methyl (CH3) group,
% but treat is as a single spin-1/2 in order to speed up simulation, we 
% need to scale the amplitude of the resulting signal by a factor of 3.
switch lower(molecule)
    case {'spin half', 'dss'}
        csVec = 0;
        JMatrix = 0;
        relAmplitudes = concentration;
        nuclei = {'1H'};
    case {'alanine', 'ala'}
        csVec   = [1.467 1.467 1.467 3.775 3.775];
        JA = 7.23; % Hz
        JMatrix = [    0     0     0    JA    JA;
                       0     0     0    JA    JA;
                       0     0     0    JA    JA;
                       JA    JA    JA    0     0;
                       JA    JA    JA    0     0];
        relAmplitudes = [1 1 1 1 1].*concentration;
        nuclei = {'1H', '1H', '1H', '1H', '1H'};
    case {'asc'} % Ascorbic acid
        %            4      5       6       6'
        csVec = [4.492      4.002   3.743   3.716];
        JMatrix = [ 0       2.07    0       0;
                    2.07    0       6       7.6;
                    0       6       0       -11.5;
                    0       7.6     -11.5   0];
        relAmplitudes = [1 1 1 1]*concentration;
        nuclei = {'1H','1H','1H','1H'};
    case 'water'
        csVec = 4.7;
        JMatrix = 0;
        relAmplitudes = 2*concentration;
        nuclei = {'1H'};
    case 'simple j'
        csVec = [-0.5 0.5]; % ppm
        J = 10; % Hz
        JMatrix = [0 J;
                   J 0]; % Hz
        relAmplitudes = [1 1].*concentration;
        nuclei = {'1H', '1H'};
    case 'naa singlet'
        csVec = 2.0;
        JMatrix = 0;
        relAmplitudes = 3*concentration;
        nuclei = {'1H'};
    case {'methanol'}
        %                   CH3              OH
        %              1      1'     1''      2
        csVec   = [  3.43    3.43    3.43    3.66];
        JMatrix = [     0       0       0     4.6;   % 1
                        0       0       0     4.6;   % 1'
                        0       0       0     4.6;   % 1''
                      4.6     4.6     4.6       0];  % 2
        relAmplitudes = [1 1 1 1]*concentration;
        nuclei = {'1H', '1H', '1H', '1H'};
    case {'eth', 'ethanol'}
        %                    CH3                  CH2        OH
        %              1      1'      1''       2     2'      3
        csVec   = [ 1.226   1.226   1.226   3.687  3.687   2.61];
        JMatrix = [     0       0       0    7.08   7.08      0;  % 1
                        0       0       0    7.08   7.08      0;  % 1' 
                        0       0       0    7.08   7.08      0;  % 1''
                     7.08    7.08    7.08       0      0      0;  % 2
                     7.08    7.08    7.08       0      0      0;  % 2'
                        0       0       0       0      0      0]; % 3
        relAmplitudes = [1 1 1 1 1 1]*concentration;
        nuclei = {'1H', '1H', '1H', '1H', '1H', '1H'};
    case {'eta', 'ethanolamine'}
        %              1      1'       2      2'
        csVec   = [ 3.818   3.818   3.147   3.147];
        JMatrix = [     0  -10.64   3.897   6.794;   % 1
                   -10.64       0   6.694   3.798;   % 1'
                    3.897   6.694       0 -11.710;   % 2
                    6.794   3.798 -11.710       0];  % 2'
        relAmplitudes = [1 1 1 1]*concentration;
        nuclei = {'1H', '1H', '1H', '1H'};
	case {'glc', 'glucose'}
		% "glc" in itself makes no reference as to whether this is the alpha
		% or the beta anomer. In the human brain, the two are in a 0.36:0.64
		% ratio; hence, both are added here (as two sub-systems) at this ratio.
        csVec   = {[5.216  3.519  3.698   3.395   3.822   3.826   3.749], [4.630  3.230  3.473   3.387   3.450   3.882   3.707]};
		JMatrix = {[    0    3.8      0       0       0       0       0;
		              3.8      0    9.6       0       0       0       0;
		                0    9.6      0     9.4       0       0       0;
		                0      0    9.4       0     9.9       0       0;
		                0      0      0     9.9       0     1.5     6.0;
		                0      0      0       0     1.5       0   -12.1;
		                0      0      0       0     6.0   -12.1       0], ...
					[   0    8.0      0       0       0       0       0;   % 1
                      8.0      0    9.1       0       0       0       0;   % 2
                        0    9.1      0     9.4       0       0       0;   % 3       
                        0      0    9.4       0     8.9       0       0;   % 4
                        0      0      0     8.9       0     1.6     5.4;   % 5
                        0      0      0       0     1.6       0   -12.3;   % 6
                        0      0      0       0     5.4   -12.3       0]};
		nuclei = {{'1H', '1H', '1H', '1H', '1H', '1H', '1H'}, {'1H', '1H', '1H', '1H', '1H', '1H', '1H'}};
		relAmplitudes = {[1 1 1 1 1 1 1]*concentration*0.36, [1 1 1 1 1 1 1]*concentration*0.64};				
    case {'glc-alpha', 'glucose alpha'}
        %              1      2      3       4       5       6      6'
        csVec   = [5.216  3.519  3.698   3.395   3.822   3.826   3.749];
        JMatrix = [    0    3.8      0       0       0       0       0;   % 1
                     3.8      0    9.6       0       0       0       0;   % 2
                       0    9.6      0     9.4       0       0       0;   % 3       
                       0      0    9.4       0     9.9       0       0;   % 4
                       0      0      0     9.9       0     1.5     6.0;   % 5
                       0      0      0       0     1.5       0   -12.1;   % 6
                       0      0      0       0     6.0   -12.1       0];  % 6'
        nuclei = {'1H', '1H', '1H', '1H', '1H', '1H', '1H'};
        relAmplitudes = [1 1 1 1 1 1 1]*concentration;
    case {'glc-beta', 'glucose-beta'}
        %              1      2      3       4       5       6      6'
        csVec   = [4.630  3.230  3.473   3.387   3.450   3.882   3.707]; 
        JMatrix = [    0    8.0      0       0       0       0       0;   % 1
                     8.0      0    9.1       0       0       0       0;   % 2
                       0    9.1      0     9.4       0       0       0;   % 3       
                       0      0    9.4       0     8.9       0       0;   % 4
                       0      0      0     8.9       0     1.6     5.4;   % 5
                       0      0      0       0     1.6       0   -12.3;   % 6
                       0      0      0       0     5.4   -12.3       0];  % 6'
        nuclei = {'1H', '1H', '1H', '1H', '1H', '1H', '1H'};
        relAmplitudes = [1 1 1 1 1 1 1]*concentration;
    case {'gly'} % Glycin
        csVec = 3.548;
        JMatrix = 0;
        relAmplitudes = 2*concentration;
        nuclei = {'1H'};
    case 'cr singlet'
        csVec = 3.0;
        JMatrix = 0;
        relAmplitudes = 3*concentration;
        nuclei = {'1H'};
    case {'cho singlet', 'ch singlet'}
        csVec = 3.2;
        JMatrix = 0;
        relAmplitudes = 9*concentration;
        nuclei = {'1H'};
    case 'cr'
        %         CH3    CH2
        csVec = {3.027, 3.913};
        relAmplitudes = {3*concentration, 2*concentration};
        JMatrix = {0, 0};
        nuclei = {{'1H'},{'1H'}};
    case 'acetate'
        csVec = 1.904;
        relAmplitudes = concentration*3;
        JMatrix = 0;
        nuclei = {'1H'};
    case 'asp'
        %              2      3     3'
        csVec   = [3.891  2.801  2.653];
        JMatrix = [    0   3.65   9.11;    % 2
                    3.65      0 -17.43;    % 3
                    9.11 -17.43      0];   % 3'
        relAmplitudes = [1 1 1]*concentration;
        nuclei = {'1H','1H','1H'};
    case 'lac'
        csVec   = [4.097  1.313  1.313  1.313];
        JMatrix = [    0   6.93   6.93   6.93;
                    6.93      0      0      0;
                    6.93      0      0      0;
                    6.93      0      0      0];
        relAmplitudes = [1 1 1 1]*concentration;
        nuclei = {'1H','1H','1H','1H'};
    case 'pcr'
        csVec = [3.029 3.930];
        relAmplitudes = concentration*[3 2];
        JMatrix = zeros(2,2);
        nuclei = {'1H','1H'};
    case {'cho', 'ch'}
        % See [1]
        %                    1      1'     2    2'
        csVec   = {3.185, [4.054 4.054 3.501 3.501]};
        JMatrix = {0, ...
                          [    0 -14.1    3.15    6.99;   % 1
                           -14.1     0    6.99    3.15;   % 1'
                            3.15  6.99       0  -14.07;   % 2
                            6.99  3.15  -14.07      0]};  % 2'
        relAmplitudes = {9*concentration, [1 1 1 1]*concentration};
        nuclei = {{'1H'},{'1H','1H','1H','1H'}};
    case {'pcho', 'pch'}
        %                      1    1'     2    2'
        csVec   = {3.209, [4.282 4.282 3.643 3.643]};
        JMatrix = {0, ...
                  [      0  -14.89    2.28    7.23;   % 1
                    -14.89       0    7.33    2.24;   % 1'
                      2.28    7.33       0  -14.19;   % 2
                      7.23    2.24  -14.19       0]};  % 2'
        relAmplitudes = {9*concentration, [1 1 1 1]*concentration};
        nuclei = {{'1H'},{'1H','1H','1H','1H'}};
    case {'gpc', 'gpcho', 'gpch'}  % Glycerophosphocholine, full
        % 
        %           ______ Choline moeity ________    _______ Glycerol moeity _____
        %          /                              \  /                             \
        csVec   = {3.212, [4.312 4.312 3.659 3.659], [3.605 3.672 3.903 3.871 3.946]};
        JMatrix = {0, ...
                   [   0 -9.32  3.10  5.90;          % 1
                   -9.32     0  5.90  3.10;          % 1'
                    3.10  5.90     0 -9.32;          % 2
                    5.90  3.10 -9.32     0],...      % 2'
                  [    0  -14.78   5.77      0       0;    % 1
                   -14.78      0   4.53      0       0;    % 1'
                    5.77    4.53      0   5.77    4.53;    % 2
                       0       0   5.77      0  -14.78;    % 3
                       0       0   4.53 -14.78       0]};  % 3'  
        relAmplitudes = {concentration*9, [1 1 1 1]*concentration, [1 1 1 1 1]*concentration};             
        nuclei = {{'1H'}, {'1H','1H','1H','1H'}, {'1H', '1H', '1H', '1H', '1H'}};
    case {'gpc-cho', 'gpc-ch'} % Glycerophosphocholine, choline moiety
        %                    1    1'     2    2'
        csVec   = {3.212, [4.312 4.312 3.659 3.659]};
        JMatrix = {0, ...
                   [   0     0  3.10  5.90;   % 1
                       0     0  5.90  3.10;   % 1'
                    3.10  5.90     0     0;   % 2
                    5.90  3.10     0     0]};  % 2'
        relAmplitudes = {concentration*9, [1 1 1 1]*concentration}; 
        nuclei = {{'1H'},{'1H','1H','1H','1H'}};
    case {'gpc-gly'}  % Glycerophosphocholine, glycine moiety
        %            1     1'    2     3     3'   
        csVec   = [3.605 3.672 3.903 3.871 3.946]; 
        JMatrix = [    0     0  5.77     0     0;   % 1
                       0     0  4.53     0     0;   % 1'
                    5.77  4.53     0  5.77  4.53;   % 2
                       0     0  5.77     0     0;   % 3
                       0     0  4.53     0     0];  % 3'
        relAmplitudes = concentration*[1 1 1 1 1];
        nuclei = {'1H', '1H', '1H', '1H', '1H'};
    case 'gaba' % See Ref. [1] for updated c.s. and J-coupling constants
        %                2       2'        3       3'        4        4'
        csVec   = [  2.284    2.284    1.889    1.889   3.0128    3.0128]; % ppm
        JMatrix = [      0  -10.744    7.755    6.173        0         0;  % 2
                   -10.744        0    7.432    7.933        0         0;  % 2'
                     7.755    7.432        0  -13.121    5.372    10.578;  % 3
                     6.173    7.933  -13.121        0    7.127     6.982;  % 3'
                         0        0    5.372    7.127        0   -12.021;  % 4
                         0        0   10.578    6.982  -12.021         0]; % 4'
        relAmplitudes = [1 1 1 1 1 1]*concentration; 
        nuclei = {'1H','1H','1H','1H','1H','1H'};
    case 'gln'
        %           3        3b       4       4b      2
        csVec   = [2.135    2.115   2.434    2.456   3.757];
        JMatrix = [    0   -14.45    6.35     6.88    5.84;   % 3
                  -14.45        0    9.16     6.88    6.53;   % 3b
                    6.35     9.16       0   -15.92       0;   % 4
                    6.88     6.88  -15.92        0       0;   % 4b
                    5.84     6.53       0        0       0];  % 2
        relAmplitudes = [1 1 1 1 1]*concentration;
        nuclei = {'1H','1H','1H','1H','1H'};
    case 'glu'
        %           3        3b       4       4b      2
        csVec   = [2.042     2.12   2.336    2.352   3.746];
        JMatrix = [    0   -14.85    6.41     8.41    7.33;   % 3
                  -14.85        0    8.48     6.88    4.65;   % 3b
                    6.41     8.48       0   -15.92       0;   % 4
                    8.41     6.88  -15.92        0       0;   % 4b
                    7.33     4.65       0        0       0];  % 2 
        relAmplitudes = [1 1 1 1 1]*concentration;
        nuclei = {'1H','1H','1H','1H','1H'};
    case 'naa'
        %            1         2      3    3'
        csVec   = {2.008, [4.382 2.6727 2.4863]};
        JMatrix = {0,...
                          [   0    3.86    9.82;    % 2
                           3.86       0  -15.59;    % 3
                           9.82  -15.59       0]};   % 3'
        relAmplitudes = {3*concentration, [1 1 1]*concentration};
        nuclei  = {{'1H'}, {'1H', '1H', '1H'}};
    case 'naag'
        % Note: the J coupling constants for the glutamate moeity are actually not
        % reported in the Govindaraju/deGraaf manuscripts. Therefore, I've used the 
        % J-coupling constants from "regular" glutamate instead.
        %
        %                   ___ Aspartate __      ___________ Glutamate __________
        %                  /                 \   /                                \
        %              1      2      3     3'      2     2       3'     4      4'
        csVec   = {2.008, [4.607  2.721  2.519], [4.128  1.881  2.049  2.190  2.180]}; 
        JMatrix = {    0, ...
                          [0      4.41   9.52;                                           % 2  \
                           4.41      0 -15.91;                                           % 3   | Aspartate
                           9.52 -15.91      0], ...                                      % 3' /
                                                 [    0 -14.85   6.41   8.41   7.33;     % 2  \
                                                 -14.85      0   8.48   6.88   4.65;     % 3   |
                                                   6.41   8.48      0 -15.92      0;     % 3'  | Glutamate
                                                   8.41   6.88 -15.92      0      0;     % 4   |
                                                   7.33   4.65      0      0      0]};   % 4' /
        relAmplitudes = {3*concentration, [1 1 1]*concentration, [1 1 1 1 1]*concentration};
        nuclei = {{'1H'},{'1H','1H','1H'},{'1H','1H','1H','1H','1H'}};
    case 'gsh'
        %          Glycine    Glutamate moeity                     Cystine moeity
        %          10         2     3       3'     4      4'       7       7'    7''     
        csVec   = {3.769,    [3.769  2.159  2.146  2.510  2.560], [4.561   2.926   2.975]};
        JMatrix = {0, ...
                  [    0   6.34   6.36      0      0;       % 2
                    6.34      0 -15.48    6.7    7.6;       % 3
                    6.36 -15.48      0    7.6    6.7;       % 3'
                       0    6.7    7.6      0 -15.92;       % 4
                       0    7.6    6.7 -15.92      0],...   % 4'
                  [    0    7.09    4.71;
                    7.09       0  -14.06;
                    4.71  -14.06       0]};     
        relAmplitudes = {2*concentration, [1 1 1 1 1]*concentration, [1 1 1]*concentration};
        nuclei = {{'1H'},{'1H','1H','1H','1H','1H'},{'1H','1H','1H'}};
    case 'gsh-glu' % Just the glutamate moeity of GSH, for speed        
        %            2     3       3'     4      4'
        csVec   = [3.769  2.159  2.146  2.510  2.560]; 
        JMatrix = [    0   6.34   6.36      0      0;    % 2
                    6.34      0 -15.48    6.7    7.6;    % 3
                    6.36 -15.48      0    7.6    6.7;    % 3'
                       0    6.7    7.6      0 -15.92;    % 4
                       0    7.6    6.7 -15.92      0];   % 4'
        relAmplitudes = [1 1 1 1 1]*concentration;
        nuclei = {'1H','1H','1H','1H', '1H'};
    case 'gsh-gly' % Just the glycine moeity of GSH, for speed
        csVec   = 3.769; 
        JMatrix = 0;
        relAmplitudes = concentration*2; 
        nuclei = {'1H'};
    case 'gsh-cys'
        %              7     7'    7''                    
        csVec   = [4.561   2.926   2.975];
        JMatrix = [    0    7.09    4.71;
                    7.09       0  -14.06;
                    4.71  -14.06       0];
        relAmplitudes = [1 1 1]*concentration;
        nuclei = {'1H', '1H', '1H'};
    case {'phosphoethanolamine', 'pe'}
        %                1       1'       2      2'
        csVec   =  [3.9765   3.9765   3.216   3.216];
        JMatrix =  [     0   -14.56   3.182   6.716;
                    -14.56        0   7.204    2.98;
                     3.216    7.204       0  -14.71;
                     6.716     2.98  -14.71       0];
        relAmplitudes = [1 1 1 1]*concentration;
        nuclei = {'1H', '1H', '1H', '1H'};
    case 'tau' % Taurine
        %            1           1'       2      2'
        csVec   =  [3.4206   3.4206  3.2459  3.2459]; 
        JMatrix =  [     0  -12.438    6.74    6.46;    % 1
                   -12.438        0    6.40    6.79;    % 1'
                      6.74     6.40       0  -12.93;    % 2
                      6.46     6.79  -12.93       0];   % 2'
        relAmplitudes = [1 1 1 1]*concentration;
        nuclei = {'1H','1H','1H','1H'};
    case {'mi', 'myo-inositol'} % Myo-inositol
        %              1      2      3      4      5      6
        csVec   = [3.522  4.054  3.522  3.614  3.269  3.614];
        JMatrix = [    0   2.89      0      0      0   10.0;   % 1
                    2.89      0   3.01      0      0      0;   % 2
                       0   3.01      0   10.0      0      0;   % 3
                       0      0   10.0      0   9.49      0;   % 4
                       0      0      0   9.49      0   9.48;   % 5
                    10.0      0      0      0   9.48      0];  % 6
        relAmplitudes = [1 1 1 1 1 1]*concentration;
        nuclei = {'1H','1H','1H','1H', '1H', '1H'};
    case {'si'} % Scyllo-Inositol
        csVec = 3.34;
        JMatrix = 0;
        relAmplitudes = 6*concentration;
        nuclei = {'1H'};
    otherwise
        error('Coupling constants not defined for molecule type %s', molecule);
end


if isempty(spins), spins = InitSpinsJ; end
if iscell(csVec)
    numSubSystems = numel(csVec);
    for idxSys=1:numSubSystems
        numNuclei = numel(csVec{idxSys});
        gmRatio{idxSys} = zeros(1,numNuclei);

        for idxNucleus=1:numNuclei
            gmRatio{idxSys}(idxNucleus) = GetGyromagneticRatio(nuclei{idxSys}{idxNucleus});
        end

        % Create equilibrium density matrix
        rho0{idxSys} = zeros(2^numNuclei, 2^numNuclei);
        for idxNucleus=1:numNuclei
            rho0{idxSys} = rho0{idxSys} + relAmplitudes{idxSys}(idxNucleus)*IzN(idxNucleus, numNuclei)/2^numNuclei;
        end

    end
else
    numNuclei = numel(csVec);
    gmRatio = zeros(1,numNuclei);

    for idxNucleus=1:numNuclei
        gmRatio(idxNucleus) = GetGyromagneticRatio(nuclei{idxNucleus});
    end

    % Create equilibrium density matrix
    rho0 = zeros(2^numNuclei, 2^numNuclei);
    for idxNucleus=1:numNuclei
        rho0 = rho0 + relAmplitudes(idxNucleus)*IzN(idxNucleus, numNuclei)/2^numNuclei;
    end
end

% Add the molecule's data
numMolecules = numel(spins.molecule);
numSpins = size(r,2);
spins.molecule(numMolecules+1).name = moleculeName; 
spins.molecule(numMolecules+1).abbrev = moleculeAbbrev; 
spins.molecule(numMolecules+1).csVec = csVec;
spins.molecule(numMolecules+1).JMatrix = JMatrix;
spins.molecule(numMolecules+1).nucleus = nuclei;
for idxSpin=1:numSpins
    spins.molecule(numMolecules+1).spin(idxSpin).r = r(:,idxSpin);
    spins.molecule(numMolecules+1).spin(idxSpin).rho = rho0;
    spins.molecule(numMolecules+1).spin(idxSpin).B1 = 1;
    spins.molecule(numMolecules+1).spin(idxSpin).B0 = 0;
    spins.molecule(numMolecules+1).spin(idxSpin).linewidth = 5;
end