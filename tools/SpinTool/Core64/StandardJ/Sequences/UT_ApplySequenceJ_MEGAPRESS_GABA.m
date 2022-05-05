% ========================================================================
% Unit Test: Pulse-Acquire for GABA using ApplySequenceJ
% ========================================================================
%
% Simulates an idealized MEGA-PRESS experiment on a simple GABA (J-coupled)
% system, using the ApplySequenceJ routine. The MEGA-PRESS sequence is 
% (following Waddell, MRI 25:1032-1038 (2007))
%
%
% FS - frequency selective editing pulse 
%
%             180     180(FS)     180      180(FS)
%      90      _        _          _          _
%       _     | |      | |        | |        | |
%      | |    | |      | |        | |        | |
% RF __| |____| |______| |________| |________| |_______[Acquire]
%       <--t1--><--t2---><---t3----><---t4----><---t5-->
%
% Following Waddell 2007, we put
%     t1 = shortest
%     t2 = TE/4
%     t3 = TE/2-t4
%     t4 = TE/4
%     t5 = TE/2-t1-t2
% such that
%     t1+t2+t3+t4+t5     = TE
%     t1-t2+t3-t4+t5     = 0   (ON)
%     t1-(t2+t3)+(t4+t5) = 0   (OFF)
%
% ========================================================================

clc
clearvars
close all

% ========================================================================
% Define spin system
% ========================================================================

csCenter = 4.7; % ppm
B0 = 2.9; % Tesla
B1 = 1; % Scaling
isSecular = 0;
linewidth = 0; % Hz

spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 
spins = SpinsJAddMolecule(spins, 'gaba');
spins = SpinsJAddMolecule(spins, 'glu');
spins = SpinsJAddMolecule(spins, 'naa');
spins.molecule(1).spin(1).linewidth = 4;
T2 = (1/(spins.linewidth + spins.molecule(1).spin(1).linewidth)/pi)*1000; % Effective T2, in ms (for ad-hoc acquisition line broadening)
% spins.molecule(2).spin(1).linewidth = 4;
% spins.molecule(3).spin(1).linewidth = 4;

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
TE = 68; % ms
t1 = 1; % Choose the "shortest possible" t1 & t2
t2 = TE/4;
t4 = TE/4;
t3 = TE/2-t4;
t5 = TE/2-t1-t2;

%                       GABA                       Glu                     NAA
%                       2.3 2.3 1.9 1.9 3.0 3.0    2.0 2.1 2.3 2.4 3.7     2.0
affectedNucleiExc = {[  1   1   1   1   1   1], [  1   1   1   1   1],    {0, [0 0 0]}};
affectedNuclei    = {[  0   0   1   1   0   0], [  1   1   1   1   0],    {1, [1 1 1]}};

seqON  = {{'hard', 90, 270, affectedNucleiExc}, {'delay', t1}, ...
          {'hard', 180, 0},                     {'delay', t2}, ...
          {'hard', 180, 0, affectedNuclei},     {'delay', t3}, ...    % Editing pulse
          {'hard', 180, 0},                     {'delay', t4}, ...
          {'hard', 180, 0, affectedNuclei},     {'delay', t5}, ...    % Editing pulse
          };
seqOFF = {{'hard', 90, 270},                    {'delay', t1}, ...
          {'hard', 180, 0},                     {'delay', t2+t3}, ...
          {'hard', 180, 0},                     {'delay', t4+t5}, ...
          };

% ========================================================================
% Apply sequence
% ========================================================================

[spinsON, ~] = ApplySequenceJ(spins, seqON);
[spinsOFF, ~] = ApplySequenceJ(spins, seqOFF);

[TTON, fidON]   = CreateTransitionTableJ(spinsON,  'isSeparateSys', true, 'AcqType', 'separate', 'Domain', 'time', 'FWHM', 6, 'TimeAxis', timeAxis);
[TTOFF, fidOFF] = CreateTransitionTableJ(spinsOFF, 'isSeparateSys', true, 'AcqType', 'separate', 'Domain', 'time', 'FWHM', 6, 'TimeAxis', timeAxis);
fprintf('Total simulation time: %.2f sec\n', toc);


%% =======================================================================
% Process and plot results
% ========================================================================

% Interleave
fidInter = {};
fidName = {};
for idxMolecule=1:numel(fidON)
    for idxSubSys=1:numel(fidON{idxMolecule})
        fidInter{end+1} = fidON{idxMolecule}{idxSubSys};
        fidInter{end+1} = fidOFF{idxMolecule}{idxSubSys};
        fidName{end+1} = sprintf('%s Sys %d (ON)', spins.molecule(idxMolecule).abbrev, idxSubSys);
        fidName{end+1} = sprintf('%s Sys %d (OFF)', spins.molecule(idxMolecule).abbrev, idxSubSys);
    end
end

spec2=SPT1DGUI('B0', B0, 'FID', fidInter, 'title', fidName, 'dwellTime', 0.001/SW);
