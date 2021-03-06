% ------------------------------------------------------------------------
% Test Script for: ApplyPulsePTX.mex
% ------------------------------------------------------------------------
% The ApplyPulsePTX function simulates the response of a spin system to
% multiple coils. We test case a few scenarios:
%
% 1. Ideal single coil transmission
%    A sinc excitation pulse is played out on a single totally homogeneous
%    coil, and a 1D sample along z.
%
% ------------------------------------------------------------------------

clear all
close all
clc

% ------------------------------------------------------------------------
% Test Case I: Ideal single coil excitation
% ------------------------------------------------------------------------

% Initialize spins: a one dimensional homogeneous distribution of thermally
% equilibrated spins (along +z), with "infinite" T1 and T2
numSpins = 80;
sampleSize = 10; % mm
numSteps = 128;
spins = InitSpinsRelax(0, numSpins, sampleSize, [0; 0; 1], 1e6, 1e6, 1);

% Create a simple sinc excitation pulse and scale its amplitude manually
% for 90-deg
sincPulse = PulseCreateSinc(4, 4, 128, 90);
sincPulse.RFamp = sincPulse.RFamp*0.53;

% Scale gradient 
sincPulse.Gz = ones(1,numSteps)*0.5;  % kHz/mm (multiply by GetGyromagneticRatio('1h')/1000 to get mT/m)

% Create an "array" of pulses for the PTX simulation. Here we only have
% one coil, so the array is kind-of redundant, but we do it for
% completeness nonetheless.
pulse(1) = sincPulse;

% Calculate the response using the single-coil simulation function
spinsOut1 = ApplyPulseRelax(spins, sincPulse);
spinsOut2 = ApplyPulsePTX(spins, pulse);

% Extract output
for k=1:numSpins
    Mz1(k) = spinsOut1(k).M(3);
    Mz2(k) = spinsOut2(k).M(3);
end

figure 
zAxis = linspace(-5, 5, numSpins);
plot(zAxis, Mz1)
hold
plot(zAxis, Mz2, 'r');
title({'M_z(z): "regular" (blue) vs. PTX Bloch solvers (red)','for a perfectly homogeneous coil (should be the same!)'});



% ------------------------------------------------------------------------
% Test Case II: Single coil excitation with an inhomogeneous B1 field
% ------------------------------------------------------------------------

% Spins structure for the multiple-coil simulation
spinsPTX = InitSpinsRelax(0, numSpins, sampleSize, [0; 0; 1], 1e6, 1e6, 1);

% Create B1 inhomogeneous profile
B1Profile = 1-(zAxis+2).^2/36;
for k=1:round(numSpins/2)
    B1Profile(k) = 1;
end
B1Profile = B1Profile.*(B1Profile>0); % Make sure no negative values
for k=1:numSpins
    spinsPTX(k).B1 = B1Profile(k);
end


% Calculate the response using the single-coil simulation function. Once
% again we use the 90-pulse we used in the previous test case.
spinsOut1 = ApplyPulseRelax(spinsPTX, sincPulse);
spinsOut2 = ApplyPulsePTX(spinsPTX, pulse);


% Extract output
for k=1:numSpins
    Mz1(k) = spinsOut1(k).M(3);
    Mz2(k) = spinsOut2(k).M(3);
end

figure
plot(zAxis, Mz1,'b')
hold
plot(zAxis, Mz2, 'r');
plot(zAxis, B1Profile, 'g');

title({'M_z(z): "regular" (blue) vs. PTX Bloch solvers (red) for an inhomogeneous single coil','Green: B1 profile (where 1 means homogeneous)'});


% ------------------------------------------------------------------------
% Test Case III: Two coils with inhomogeneous B1 field distributions
% ------------------------------------------------------------------------

% Spins structure for the multiple-coil simulation
spinsPTX = InitSpinsRelax(0, numSpins, sampleSize, [0; 0; 1], 1e6, 1e6, 1);

% Create B1 inhomogeneous profile
B1Profile1 = 1-(zAxis+2).^2/36;
for k=1:round(numSpins/2)
    B1Profile1(k) = 1;
end
B1Profile1 = B1Profile1.*(B1Profile>0); % Make sure no negative values

% The second B1 profile is the mirror image of the first
B1Profile2 = B1Profile1(end:-1:1);

for k=1:numSpins
    spinsPTX(k).B1 = [B1Profile1(k) B1Profile2(k)];
end

% Create the two channel pulse: actually, just transmit both pulse on
% both channels
pulsePTX(1) = sincPulse;
pulsePTX(2) = sincPulse;

% Calculate the response using the single-coil simulation function. Once
% again we use the 90-pulse we used in the previous test case.
spinsOut1 = ApplyPulseRelax(spinsPTX, sincPulse);
spinsOut2 = ApplyPulsePTX(spinsPTX, pulsePTX);


% Extract output
for k=1:numSpins
    Mz1(k) = spinsOut1(k).M(3);
    Mz2(k) = spinsOut2(k).M(3);
end

figure
plot(zAxis, Mz1,'b')
hold
plot(zAxis, Mz2, 'r');
plot(zAxis, B1Profile1, 'g');
plot(zAxis, B1Profile2, 'g.');

title({'M_z(z): "regular" (blue) vs. PTX Bloch solvers (red) for an inhomogeneous single coil','Green: B1 profiles of both coils'});
