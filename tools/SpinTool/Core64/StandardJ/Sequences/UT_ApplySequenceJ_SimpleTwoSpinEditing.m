% ========================================================================
% Unit Test: Spectral editing for a simple two-spin IS system
% ========================================================================
%
% Simulates an idealized editing sequence for an IS system, of the form
%
%     90(-y)      180(x),S    180(x)      180(x),S
%                 On/Off                  On/Off
%                   _           _           _
%       _          | |         | |         | |
%      | |         | |         | |         | |
% RF __| |_________| |_________| |_________| |__________[Acquire]
%        <------- TE/2 ---------><-------- TE/2 ------->
%
% On the OFF scan, we have a simple spin echo, and the I spin evolves as:
%
%    90(-y)      TE=1/J
% Iz ------> Ix --------> Ix*cos(pi*J*TE) + 2*IyIz*sin(pi*J*TE) = -Ix
% 
% While on the ON scan, the I gets completely refocused:
%
%    90(-y)        TE
% Iz ------> Ix --------> Ix
%  
% On the other hand, a non-coupled system will be refocused on both counts,
% so subtracting OFF - ON will cancel out the singlet but enhance the I
% resonance two-fold.
% ========================================================================

clear all
close all
clc

% ========================================================================
% Define spin system
% ========================================================================

csCenter = 0; % ppm
B0 = 3; % Tesla
B1 = 1; % Scaling
isSecular = 1;
linewidth = 0; % Hz

spins = InitSpinsJ(csCenter, B0, isSecular, linewidth, B1); 
spins = SpinsJAddMolecule(spins, 'Simple J');
spins = SpinsJAddMolecule(spins, 'spin half');
spins.molecule(1).spin(1).linewidth = 4;
spins.molecule(2).spin(1).linewidth = 4;
J = max(spins.molecule(1).JMatrix(:)); % Hz

% ========================================================================
% Define editing sequence
% ========================================================================

Gx = 0;
Gy = 0;
Gz = 0;
SW = 1.0; % kHz
numAcqPoints = 1000; 
TE = 1/J*1000; % ms

affectedNuclei = {[0 1], [1]};

seqON  = {{'hard', 90, 270},                {'delay', TE/4}, ...
          {'hard', 180, 0},                 {'delay', TE/4}, ...
          {'hard', 180, 0, affectedNuclei}, {'delay', TE/4}, ...  
          {'hard', 180, 0},                 {'delay', TE/4}, ...
          {'acquire', numAcqPoints, SW, Gx, Gy, Gz}};
seqOFF = {{'hard', 90, 270},                {'delay', TE/2}, ...
          {'hard', 180, 0},                 {'delay', TE/2}, ...
          {'acquire', numAcqPoints, SW, Gx, Gy, Gz}};

% ========================================================================
% Apply sequence
% ========================================================================

[spinsON, fidON] = ApplySequenceJ(spins, seqON);
[spinsOFF, fidOFF] = ApplySequenceJ(spins, seqOFF);

% ========================================================================
% Process and plot results
% ========================================================================

spec=SPT1DSpectrum('B0', 3.0, 'FID', {fidON{1}, fidOFF{1}}, 'dwellTime', 0.001/SW);
