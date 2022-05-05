% ================================================================
% Example Script: How to use 3D spin structures
% ================================================================
%
% 3D spin structures are used to describe a spin distribution on
% a rectilinear grid. All spins share the same "properties", meaning
% T1, T2, chemical shift, equilibrium magnetization (M0), etc ...
% This allows for considerable speedups when simulation 3D spin
% ensembles of non-coupled spin-1/2s.
%
% NOTE: This structure should ONLY BE USED with constant-gradient
% pulses!
%
% ================================================================

clear all
close all
clc


% InitSpins3D is used to initialize the spin structure. It accepts
% parameter name - value pairs.
numSpins = [81 81 81];
spins = InitSpins3D('numSpins', numSpins, ...     % Along x, y, z axes
                    'sampleSize', [40 40 60], ...   % In mm
                    'T1', 1500);

% Create a 200 Hz 90-deg Gaussian excitation pulse. Add a constant
% gradient along z-axis.
pulse = PulseCreateGaussian(0, 0.2, 2, 90);       
pulse.Gz = pulse.Gz + 0.01; % kHz/mm
seq = {pulse};

% 3D spins are propagated using the ApplySequence3D command.
tic
fprintf('Running 3D simulation ... \n');
spinsOut = ApplySequence3D(spins, seq);
fprintf('Simulation ended. Total time needed: %.2f sec\n', toc);

% For the sake of comparison, let's run the same thing using the
% full spin structure
spinsFull = InitPhantomCube([40 40 60], [0 0 0], numSpins, 1500, 200, 1, 0, 1);
fprintf('Running conventional simulation ... \n');
spinsFullOut = ApplySequence(spinsFull, seq);
fprintf('Simulation ended. Total time needed: %.2f sec\n', toc);

% To extract the simulated magnetization from the 3D spin structure,
% we need to access the magnteization directly. It is stored in the
% 4D array spinsOut.M (x*y*z*3)
Mz = spinsOut.M(:,:,:,3);
MzSlice = squeeze(Mz(round(numSpins(1)/2),:,:)); % Take a slice from yz plane
figure
subplot(1,2,1)
imagesc(spinsOut.zVec, spinsOut.yVec, MzSlice)
title('M_z (yz slice)');
xlabel('z (mm)');
ylabel('y (mm)');
subplot(1,2,2)
plot(spinsOut.zVec, squeeze(MzSlice(round(numSpins(2)/2), :)));
title('Plot along y-const line');
xlabel('z (mm)');
ylabel('M_z');
ylim([-1 1]);
