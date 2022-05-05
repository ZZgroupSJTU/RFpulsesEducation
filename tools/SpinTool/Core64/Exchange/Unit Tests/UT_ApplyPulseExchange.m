% ========================================================================
% Creates two exchanging peaks ("solvent" and "solute") and plots the
% z-spectrum, i.e. the size of the water peak as a function of the 
% irradiation frequency, for the given irradiation pattern.
% ========================================================================

clear all
close all
clc

isTest1 = 0;
isTest2 = 0;

fprintf('This simulation will test the function ApplyPulseExchange.m.\n');
fprintf('\n');
x=input('Hit Enter to continue ...');
fprintf('\n');

% ========================================================================
% Test 1: Two precessing uncoupled spins, irradiation @ 0 kHz
% ========================================================================

if isTest1

    fprintf('=============================================================\n');
    fprintf('Unit Test 1\n');
    fprintf('=============================================================\n');
    fprintf('First, we create two non-coupled spins and apply a 1ms 90(-y)\n');
    fprintf('pulse and see that it rotates our magnetization from the z-axis \n');
    fprintf('to the x-axis.\n');
    fprintf('\n');
    x=input('Hit Enter to continue ...');
    fprintf('\n');

    % Create spins
    spinOffset = 0.3; % kHz
    chemShiftVec = [0 spinOffset];  % kHz
    k = 0; % Hz
    exchangeMatrix = [0 k; -k 0];
    initialMag = [0; 0; 1; 0; 0; 1];
    T1 = 1e9;
    T2 = 1e9; 
    eqMag = [1 1];
    spins = InitSpinsExchange(chemShiftVec, exchangeMatrix, initialMag, T1, T2, eqMag);

    pulseOffsetPPM = 0;
    pulseDuration = 1; % ms
    pulsePhase = 3*pi/2; % radians
    flipAngle = 90; % degrees
    peakB1 = flipAngle/360*1/pulseDuration; % kHz
    numSteps = 2;
    pulseRect = PulseCreateConst(pulseDuration, numSteps, peakB1, pulsePhase);
    pulseRect.offset = 0;

    fprintf('Input: \n');
    fprintf('   Spins: Chemical shifts: [0 0.3] kHz \n');
    fprintf('          T1=T2=infinity \n');
    fprintf('   Pulse: 1 ms at gammabar*B1=0.25 kHz \n');
    fprintf('          %d time steps \n', numSteps);
    fprintf('          Offset: 0 kHz \n');
    fprintf('Predicted final magnetization: \n');
    fprintf('   Spins: Chemical shifts: [0 0.3] kHz \n');
    fprintf('   Pulse: 1 ms at gammabar*B1=0.25 kHz \n');

    fprintf('Initial magnetization: \n');
    fprintf('   Spin #1: (%.5f, %.5f, %.5f) \n', spins(1).M(1), spins(1).M(2), spins(1).M(3));
    fprintf('   Spin #2: (%.5f, %.5f, %.5f) \n', spins(1).M(4), spins(1).M(5), spins(1).M(6));
    fprintf('\n');

    spins = ApplyPulseExchange(spins, pulseRect);
    fprintf('Final magnetization: \n');
    fprintf('   Spin #1: (%.5f, %.5f, %.5f), Magnitude: %.5f \n', spins(1).M(1), spins(1).M(2), spins(1).M(3), norm(spins(1).M(1:3)) );
    fprintf('   Spin #2: (%.5f, %.5f, %.5f), Magnitude: %.5f \n', spins(1).M(4), spins(1).M(5), spins(1).M(6), norm(spins(1).M(4:6)));
    fprintf('\n');
    fprintf('Predicted magnetization (based on rotation matrices): \n');
    M1 = RotMat([0 -peakB1 0], peakB1*pulseDuration*2*pi)*[0;0;1];
    M2 = RotMat([0 -peakB1 spinOffset], sqrt(peakB1^2 + spinOffset^2)*pulseDuration*2*pi)*[0;0;1];
    fprintf('   Spin #1: (%.5f, %.5f, %.5f), Magnitude: %.5f \n',  M1(1), M1(2), M1(3), norm(M1));
    fprintf('   Spin #2: (%.5f, %.5f, %.5f), Magnitude: %.5f \n',  M2(1), M2(2), M2(3), norm(M2));
    fprintf('\n');
    x=input('Hit Enter to continue ...');
    fprintf('\n');
end

% ========================================================================
% Test 2: Two precessing uncoupled spins, irradiation @ 0.3 kHz
% ========================================================================

if isTest2
    fprintf('=============================================================\n');
    fprintf('Unit Test 2\n');
    fprintf('=============================================================\n');
    fprintf('Now we repeat the previous example, setting the RF offset to \n');
    fprintf('match that of the second spin. Exchange constant is 0.\n'); 
    fprintf('\n');
    x=input('Hit Enter to continue ...');
    fprintf('\n');

    % Create spins
    spinOffset = 0.3; % kHz
    chemShiftVec = [0 spinOffset];  % kHz
    k = 0; % Hz
    exchangeMatrix = [0 k; -k 0];
    initialMag = [0; 0; 1; 0; 0; 1];
    T1 = 1e9;
    T2 = 1e9; 
    eqMag = [1 10];
    spins = InitSpinsExchange(chemShiftVec, exchangeMatrix, initialMag, T1, T2, eqMag);

    pulseDuration = 1; % ms
    pulsePhase = 3*pi/2; % radians
    flipAngle = 90; % degrees
    peakB1 = flipAngle/360*1/pulseDuration; % kHz
    numSteps = 8;
    pulseRect = PulseCreateConst(pulseDuration, numSteps, peakB1, pulsePhase);
    pulseRect.offset = spinOffset; % kHz

    fprintf('Input: \n');
    fprintf('   Spins: Chemical shifts: [0 0.3] kHz \n');
    fprintf('          T1=T2=infinity \n');
    fprintf('   Pulse: 1 ms at gammabar*B1=0.25 kHz \n');
    fprintf('          %d time steps \n', numSteps);
    fprintf('          Offset: 0 kHz \n');
    fprintf('Predicted final magnetization: \n');
    fprintf('   Spins: Chemical shifts: [0 0.3] kHz \n');
    fprintf('   Pulse: 1 ms at gammabar*B1=0.25 kHz \n');

    fprintf('Initial magnetization: \n');
    fprintf('   Spin #1: (%.5f, %.5f, %.5f) \n', spins(1).M(1), spins(1).M(2), spins(1).M(3));
    fprintf('   Spin #2: (%.5f, %.5f, %.5f) \n', spins(1).M(4), spins(1).M(5), spins(1).M(6));
    fprintf('\n');

    spins = ApplyPulseExchange(spins, pulseRect);
    fprintf('Final magnetization: \n');
    fprintf('   Spin #1: (%.5f, %.5f, %.5f), Magnitude: %.5f \n', spins(1).M(1), spins(1).M(2), spins(1).M(3), norm(spins(1).M(1:3)) );
    fprintf('   Spin #2: (%.5f, %.5f, %.5f), Magnitude: %.5f \n', spins(1).M(4), spins(1).M(5), spins(1).M(6), norm(spins(1).M(4:6)));
    fprintf('\n');
    fprintf('Predicted magnetization (based on rotation matrices): \n');
    T = RotMat([0 0 1], -pulseRect.offset*pulseDuration*2*pi); % Transformation matrix from pulse's frame to "original frame"
    R1 = RotMat([0 -peakB1 0-pulseRect.offset], sqrt(peakB1^2+pulseRect.offset^2)*pulseDuration*2*pi); % Rotation in "pulse rotating frame" for first spin (now at -0.3 kHz)
    R2 = RotMat([0 -peakB1 spinOffset-pulseRect.offset], peakB1*pulseDuration*2*pi); % Rotation in "pulse rotating frame" for first spin (now at 0 kHz)
    M1 = T*R1*[0;0;1];
    M2 = T*R2*[0;0;1];
    fprintf('   Spin #1: (%.5f, %.5f, %.5f), Magnitude: %.5f \n',  M1(1), M1(2), M1(3), norm(M1));
    fprintf('   Spin #2: (%.5f, %.5f, %.5f), Magnitude: %.5f \n',  M2(1), M2(2), M2(3), norm(M2));
end
