clear all 
close all
clc

% Let's simulate a simple 90-acquire experiment

% Step 1: Define the spin structure
chemShiftVec = [0.05 0.25];  % One chemical shift, at w=0 kHz
numSpinsPerShift = 1;  % Only need one spin since you don't have imaging
sampleSize = 100; % Irrelevant, only one spin
initMag = [0; 0; 1];
T1 = 1000; % ms
T2 = 50; % ms
eqMag = 1; % M0
sampleOffset = 0; % mm. Irrelevant, only one spin.
B1 = 1.0; % This has to do with transmission inhomogeneity. Set this to 1 (disregard for now).
spatialAxis = 'z';  % Doesn't matter. Can also be 'x' or 'y'.
spins = InitSpinsRelax(chemShiftVec, numSpinsPerShift, sampleSize, initMag, T1, T2, eqMag, sampleOffset, B1, spatialAxis);

% Step 2: Define the sequence
tiltAngle = 90;
pulsePhase = 0;

numAcqPts = 512;
SW = 1; % kHz 
% The spectral width SW determines the dwell time, which is the time between
% successive points. SW = 1 kHz means your dwell time is 1/SW = 1 ms.
% After a Fourier transform, what you will see is that the range of frequencies
% is indeed SW (that is, the dwell time determines the range of frequencies
% you can observe).
seq = {{'hard', tiltAngle, pulsePhase}, {'acquire', numAcqPts, SW, 0, 0, 0}};

% Step 3: Apply the sequence
[spinsOut,fidCellArray]=ApplySequence(spins, seq);

% Step 4: Let's plot the FID
figure
subplot(2,1,1);
timeAxis = [0:(1/SW):(numAcqPts-1)*(1/SW)]; 
plot(timeAxis, real(fidCellArray{1}));
hold on
plot(timeAxis, imag(fidCellArray{1}), 'r');
legend({'M_x', 'M_y'});
xlabel('Time (ms)');
ylabel('Magnetization (a.u.)');
title('FID');

% Step 5: Plot the spectrum
spec = fftshift(fft(fidCellArray{1}));
dv = SW/numAcqPts;
freqAxis = [-SW/2:dv:(SW/2-dv)];
subplot(2,1,2);
plot(freqAxis, real(spec));
hold on
plot(freqAxis, imag(spec), 'r');
xlabel('Freq. (kHz)');
ylabel('Spectrum (a.u.)');
legend({'Real', 'Imag'});
title('Spectrum');

% TODO:
% 1. In the basic 90-acquire experiment, play around with the following:
%    1a. What is the effect of the tilt angle on the signal?
%    1b. What is the effect of taking cs>SW? (i.e. "observe" aliasing/wrapping)
%    1c. How does the pulsePhase affect the signal/spectrum?
%    1d. How does T2 affect your spectrum? 
%    1e. Why doesn't T1 affect your spectrum?
%    1f. How does the number of points (or, alternatively, SW) affect your spectrum?
% 2. Implement an inversion recovery experiment:
%    2a. Play around with TI and see that you get an inversion recovery
%        (that is, that the signal intensity diminishes and then increases)
%    2b. How should you choose TI to make your signal vanish?




return
%%
clear all
close all
clc
% Variable Name    Units    Description
% frequencies      Hz       Vector containing frequencies
% intensities      a.u.     Intensity of each peak. A vector or number [1].
% T2               ms       Transverse relaxation. Either a vector or a 
%                           number [1]. 
% acqTime          ms       Total acquisition time
% numPoints        -        Number of acquisition points
% noiseSTD         a.u.     STD of (Gaussian!) noise (set to 0 for none)
% ZF               -        Amount of zero filling (set to 0 for no zero
%                           filling), in multiples of numPoints
SW = 1; % kHz
frequencies = [-250 -125 -80 10 50 55 66 100];
intensities = [   1    1   1  1  1  1  1   1];
T2 = 12;
noiseSD = 0.01;
ZF = 0;
numPoints = 512;
dwellTime = 1/SW;
acqTime = dwellTime*numPoints; % ms
DCOffsetSpec = 0;
[tt, fid, vv, spec, ~, ~, basisSpec] = GenerateDummyFID(frequencies, intensities, T2, acqTime, numPoints, noiseSD, ZF, DCOffsetSpec);
BW = 1/dwellTime; % kHz

[frequenciesSVD, dampingsSVD, basisSVD, ahatSVD] = lpsvd3(fid.', SW);
fidSVD = basisSVD*ahatSVD;
specSVD = fftshift(fft(fidSVD));

% Now let's reconstruct the spectrum using only the numFiltCoeffs largest 
% components in the SVD. To do so, we first rearrange the coefficients
% in descending order according to magnitude, and then take the first
% numFiltCoeffs ones.
numFiltCoeffs = min(20, numel(ahatSVD));
[~, indices] = sort(abs(ahatSVD), 'descend');
ahatSVDSorted = ahatSVD(indices); % Don't use the output of sort - it spits out the absolute ahat values!
basisSVDSorted = basisSVD(:, indices);
fidSVDFilt = basisSVDSorted(:,1:numFiltCoeffs)*ahatSVDSorted(1:numFiltCoeffs);
specSVDFilt = fftshift(fft(fidSVDFilt));

% Finally, let's create a spectrum where we remove certain components. 
frequenciesSVDSorted = frequenciesSVD(indices)*1000; % Convert to Hz
removeRange = [30 70]; % Hz
indicesToRemove = intersect( find( frequenciesSVDSorted>= removeRange(1)), find( frequenciesSVDSorted<= removeRange(2)));
fidSVDSortedFilt = basisSVDSorted(:, setdiff(1:numFiltCoeffs, indicesToRemove))*ahatSVDSorted(setdiff(1:numFiltCoeffs, indicesToRemove));
specSVDSortedFilt = fftshift(fft(fidSVDSortedFilt));

%%
figure
ax1=subplot(4,1,1);
plot(vv, real(spec))
hold on
plot(vv, real(specSVD), 'r--');
plot(vv, real(specSVDFilt), 'm.-');
plot(vv, real(specSVDSortedFilt), 'k.-');
title('Spectrum & Reconstructed SVD Spectrum');
legend({'Spectrum', 'SVD Spectrum', '"Denoised" SVD Spectrum"', 'Filtered denoised SVD'});
xlabel('Hz');
ylabel('a.u.');

ax2=subplot(4,1,2);
basisSVDSpec = fftshift(fft(basisSVD, [], 1), 1);
plot(vv, real(basisSVDSpec.*repmat(ahatSVD.', [numPoints, 1])));
title('Individual basis functions: SVD');
xlabel('Hz');
ylabel('a.u.');

ax3=subplot(4,1,3);
basisSVDSpecSorted = fftshift(fft(basisSVDSorted(:, 1:numFiltCoeffs), [], 1), 1);
plot(vv, real(basisSVDSpecSorted.*repmat(ahatSVDSorted(1:numFiltCoeffs).', [numPoints, 1])));
title(sprintf('Individual basis functions: SVD (%d largest coeffs)', numFiltCoeffs));
xlabel('Hz');
ylabel('a.u.');

ax4=subplot(4,1,4);
plot(vv, real(basisSpec));
title('Individual basis functions: True');
xlabel('Hz');
ylabel('a.u.');

linkaxes([ax1, ax2, ax3, ax4]);