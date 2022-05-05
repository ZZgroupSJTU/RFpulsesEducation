% ========================================================================
% Unit Test: TE-Averaged PRESS for Glutamate Detection
% ========================================================================
%
% Simulates an idealized TE-Averaged PRESS for the detection of Glutamate
% following Hurd, Magn Reson Med 51:435-440 (2004)
%
% FS - frequency selective editing pulse 
%
%             180          180      
%      90      _            _         
%       _     | |          | |        
%      | |    | |          | |       
% RF __| |____| |__________| |_______[Acquire]
%       <-----TE1----><----TE2------>
%
% ========================================================================

% clear all
% close all
clc

% ========================================================================
% Define spin system
% ========================================================================

csCenter = 4.7; % ppm
B0 = 2.9; % Tesla
B1 = 1; % Scaling
isSecular = 0;
linewidth = 0; % Hz

spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 
spins = SpinsJAddMolecule(spins, 'glu');
% spins = SpinsJAddMolecule(spins, 'gln');
spins = SpinsJAddMolecule(spins, 'naa singlet');
spins.molecule(1).spin(1).linewidth = 4;
T2 = (1/(spins.linewidth + spins.molecule(1).spin(1).linewidth)/pi)*1000; % Effective T2, in ms (for ad-hoc acquisition line broadening)
spins.molecule(2).spin(1).linewidth = 4;

% ========================================================================
% Define editing sequence:
% ========================================================================

Gx = 0;
Gy = 0;
Gz = 0;
SW = 1.2; % kHz
numAcqPoints = 1000; 
dt = 1/SW;
timeAxis = [0:dt:(numAcqPoints-1)*dt];
TE1 = 15;
TE2 = 15:10:110;

affectedNucleiExc = {[1 1 1 1 1 1]};

% ========================================================================
% Apply sequence
% ========================================================================

fprintf('Simulating TE2 (ms) = ');
for idxTE2=1:numel(TE2)
    fprintf('%.1f   ', TE2(idxTE2));
    seq = {{'hard', 90, 270},                    {'delay', TE1/2}, ...
           {'hard', 180, 0},                     {'delay', TE1/2+TE2(idxTE2)/2}, ...
           {'hard', 180, 0},                     {'delay', TE2(idxTE2)/2}, ...
           };
    [spinsOut, ~] = ApplySequenceJ(spins, seq);
    TT = CreateTransitionTableJ(spinsOut);
    fid{idxTE2} = zeros(1, numAcqPoints);    
    numLines = size(TT,1);
    for idx=1:numLines
        fid{idxTE2} = fid{idxTE2} + TT(idx,2)*exp(-timeAxis/T2).*exp(-1i*2*pi*TT(idx,1)*timeAxis);
    end
end
fprintf('Done!\n');


% ========================================================================
% Process and plot results
% ========================================================================

% All you have to do now is add up all the different TEs!
spec2=SPT1DGUI('B0', B0, 'FID', fid, 'dwellTime', 0.001/SW);
