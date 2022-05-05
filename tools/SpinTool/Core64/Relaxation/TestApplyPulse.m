% ========================================================================
% 1D Gradient Echo, by Assaf Tal (Jul 2011)
% ------------------------------------------------------------------------
%
% This script uses the ApplyPulse function to simulate and then post
% -process a simple 1D GRE experiment. The comments provide the necessary 
% background information.
%
% ========================================================================

clear all
close all
clc

% ------------------------------------------------------------------------
% Create sample
% ------------------------------------------------------------------------

% There are two basic structures used in every simulation: spins and
% pulses. A spin structure is an array, with each entry representing a
% different spin. Each spin has several pre-defined properties: position,
% chemical shift, T1, T2, etc. A spin structure is your "phantom", so to
% speak. 
%
% For this simulation we'll create a simple "boxcar" phantom in 1D.
%
% Note: The fact that it ISN'T a cartesian 3D matrix allows you to
% efficiently create non-cartesian and/or sparse phantoms (i.e. not waste 
% time simulating empty space).
numSpins = 10000;
M0 = 1; % Equilibrium magnetization, in arbitrary units

% We'll take our sample to be along the z-axis (with x=y=0). We now
% calculate the z ordinate of each spin, in mm. This is done out of
% convenience at this point so we won't have to do it when populating the 
% spin structure.
sampleSize = 100; % In mm
spinPositions = linspace(-sampleSize/2, sampleSize/2, numSpins);


for idxSpin=1:numSpins
    % We now proceed to define, for the idxSpin spin, all relevant
    % characteristics:
    % Define the spin's position in mm. 
    spins(idxSpin).r = [0; 0; spinPositions(idxSpin)];
    % Define the magnetization vector of each spin. We assume we start from
    % thermal equilibrium, which is the same for all spins.
    spins(idxSpin).M = [0; 0; M0];
    % Define the chemical shift in kHz
    spins(idxSpin).cs = 0;
    % Define T1 and T2 in ms
    spins(idxSpin).T1 = 1000;
    spins(idxSpin).T2 = 1000;
    % Determine equilibrium value of magnetization
    spins(idxSpin).M0 = 1;
    % Determine B0, B1-transmit and B1-receive Inhomogeneities
    spins(idxSpin).B0 = 0; % An offset in kHz.
    spins(idxSpin).B1 = 1; % Scaling factor during transmit
    spins(idxSpin).RS = 1; % Scaling factor during receive
end



% ------------------------------------------------------------------------
% Create excitation pulse
% ------------------------------------------------------------------------
%
% Next we're going to define the pulse the will excite the sample.
% This means using the "Pulse" structure. This basically holds the RF &
% phase of the pulse in the form of arrays, as well as any accompanying
% gradients. For a 1D GRE I'd just define a hard pulse as follows:

numPulseSteps = 5000;
tp = 1; % in ms
pulseHard.tp = tp; % (in ms)
pulseHard.RFamp = 1/(4*tp) * ones(1,numPulseSteps); % (in kHz)
pulseHard.RFphase = zeros(1, numPulseSteps); % in radians
pulseHard.Gx = zeros(1, numPulseSteps);  % in kHz/mm, i.e., take the mT/m and multiply by
                   % GetGyromagneticRatio('1h')/1000 for protons
pulseHard.Gy = zeros(1, numPulseSteps);  % in kHz/mm
pulseHard.Gz = zeros(1, numPulseSteps);  % in kHz/mm

% This is a somewhat degenerate case since the hard pulse only has one
% step. In general, all fields (except tp) can be arrays. However, they all
% need to have the same number of elements.


% ------------------------------------------------------------------------
% Create acquisition "pulse"
% ------------------------------------------------------------------------
%
% This is not a pulse per-se, since RFamp = 0, but it will hold the
% gradient values during acquisition. For simplicity we assume the gradient
% doesn't have any ramp up or down times (instantaneous).

% Define acquisition parameters
numVoxels = 50; 
FOV = sampleSize*2; % in mm. We want something larger than the sample.

dK = 1/FOV; % in mm^(-1)
kMax = dK*numVoxels/2;  % We scan k-space symmetrically

acqTime = 10; % Acquisition time - short! 10 ms.
acqGradient = (kMax*2)/acqTime; % in kHz/mm
numAcqSteps = numVoxels; 

pulseAcq.tp = acqTime; % Total acquisition time
pulseAcq.RFamp = zeros(1,numAcqSteps);
pulseAcq.RFphase = zeros(1,numAcqSteps);
pulseAcq.Gx = zeros(1,numAcqSteps);
pulseAcq.Gy = zeros(1,numAcqSteps);
pulseAcq.Gz = acqGradient*ones(1,numAcqSteps);

% ------------------------------------------------------------------------
% Create rewinding gradient "pulse"
% ------------------------------------------------------------------------
%
% We need to "rewind" the gradient (go to -kMax) before we start.

pulseRewind.tp = acqTime/2; % Total acquisition time
pulseRewind.RFamp = 0;
pulseRewind.RFphase = 0;
pulseRewind.Gx = 0;
pulseRewind.Gy = 0;
pulseRewind.Gz = -acqGradient;


% ------------------------------------------------------------------------
% Finally, we're ready to simulate! This is the easy part.
% ------------------------------------------------------------------------

tic
spins = ApplyPulseRelax(spins, pulseHard);
spins = ApplyPulseRelax(spins, pulseRewind);
[spins,fid] = AcquirePulseRelax(spins, pulseAcq);
fprintf('Total time elapsed: %.2f sec \n', toc);
% That's it!

% ------------------------------------------------------------------------
% Post-process and display image
% ------------------------------------------------------------------------

y = fftshift(fft(fid));

plot(abs(y));  % I'll let you figure out how to phase correct y so you can display real(y)!
