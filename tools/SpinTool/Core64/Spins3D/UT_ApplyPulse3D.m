clc
clear all
close all
% mex ApplyPulse3D.c;

% Unit Test #1
spins = InitSpins3D('numSpins', [1 1 1]);
pulse = PulseCreateHard(1, 90, 270, 1);
spinsOut = ApplyPulse3D(spins, pulse);
fprintf('UT1: Applying a single step 90(-y) hard pulse to a single spin starting from z-axis.\n');
fprintf('Resulting magnetization = [%.3f %.3f %.3f] (should be [1 0 0])\n\n', spinsOut.M(1), spinsOut.M(2), spinsOut.M(3));

% Unit Test #2
spins = InitSpins3D('numSpins', [1 1 1]);
pulse = PulseCreateHard(1, 90, 90, 10);
spinsOut = ApplyPulse3D(spins, pulse);
fprintf('UT2: Applying a 10-step 90(y) pulse to a single spin starting from z-axis.\n');
fprintf('Resulting magnetization = [%.3f %.3f %.3f] (should be [-1 0 0])\n\n', spinsOut.M(1), spinsOut.M(2), spinsOut.M(3));


% Unit Test #3
spins = InitSpins3D('numSpins', [3 5 7]);
pulse = PulseCreateHard(1, 90, 90, 10);
spinsOut = ApplyPulse3D(spins, pulse);
fprintf('UT3: Applying a 10-step 90(y) pulse to a 3D array of spins starting from z-axis.\n\n');
figure
imagesc(squeeze(spinsOut.M(:,:,3,1)));
title({'UT3: M_x in xy slice after pulse','(Should be -1 everywhere)'});

% Unit Test #4
fprintf('UT4: Applying an SLR90 pulse along x-axis to a 3D array of spins starting from z-axis.\n');
numSpins = [30 31 32];
sampleSize = [30 30 30];
spins = InitSpins3D('numSpins', numSpins, 'sampleSize', sampleSize);
pulse = PulseCreateSLR90(100, 1, 2, 'linearr6', 0, 10, 'x');
tic
spinsOut = ApplyPulse3D(spins, pulse);
fprintf('Time: %.3f sec \n\n', toc);
figure
subplot(1,2,1)
imagesc(spinsOut.yVec, spinsOut.xVec, squeeze(spinsOut.M(:,:,15,1)));
xlabel('y (mm)')
ylabel('x (mm)')
title({'UT4: M_x in xy slice after pulse','(Should excite 10 mm slice)'});
subplot(1,2,2)
plot(spinsOut.xVec, squeeze(spinsOut.M(:,16,15,1))); 
response = PlotPulseFreqResponse(pulse, [0;0;1], -sampleSize(1)/2, sampleSize(1)/2, numSpins(1), 'mx', 'x');
hold on
plot(spinsOut.xVec, response{1}, 'r--')
xlabel('x (mm)');
ylabel('Mx (for given y, z)');
legend({'ApplyPulse3D', 'ApplyPulseRelax'});

% Unit Test #5
fprintf('UT5: Applying an SLR90 pulse along x-axis to a 3D array of spins starting from x-axis.\n');
numSpins = [30 31 32];
sampleSize = [30 30 30];
spins = InitSpins3D('numSpins', numSpins, 'sampleSize', sampleSize, 'initMag', [1;0;0]);
pulse = PulseCreateSLR90(100, 1, 2, 'linearr6', 0, 10, 'x');
tic
spinsOut = ApplyPulse3D(spins, pulse);
fprintf('Time: %.3f sec \n', toc);
figure
subplot(1,2,1)
imagesc(spinsOut.yVec, spinsOut.xVec, squeeze(spinsOut.M(:,:,15,1)));
xlabel('y (mm)')
ylabel('x (mm)')
title({'UT5: M_x in xy slice after pulse','(Should excite 10 mm slice)'});
subplot(1,2,2)
plot(spinsOut.xVec, squeeze(spinsOut.M(:,16,15,1))); 
response = PlotPulseFreqResponse(pulse, [1;0;0], -sampleSize(1)/2, sampleSize(1)/2, numSpins(1), 'mx', 'x');
hold on
plot(spinsOut.xVec, response{1}, 'r--')
xlabel('x (mm)');
ylabel('Mx (for given y, z)');
legend({'ApplyPulse3D', 'ApplyPulseRelax'});

%% Unit Test #6
fprintf('\nUT6: Testing relaxation. Start from single spin, M=[1 0 0] and measure after a time 1 sec.\n');
spins = InitSpins3D('numSpins', [1 1 1], 'initMag', [1 0 0]);
T = 200;
pulse = PulseCreateZero(T, 1); % Zero pulse, 1 sec long, with 1 step
spinsOut = ApplyPulse3D(spins, pulse);
fprintf('Resulting magnetization = [%.3f %.3f %.3f] (should be [%.3f %.3f %.3f])\n\n', spinsOut.M(1), spinsOut.M(2), spinsOut.M(3), 1*exp(-T/spinsOut.T2), 0, (1-exp(-T/spinsOut.T1)));

