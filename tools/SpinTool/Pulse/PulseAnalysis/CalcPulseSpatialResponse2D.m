function [vecPosX, vecPosY, Mz, Mx, My] = CalcPulseSpatialResponse2D(pulse,...
                                            zCoordinate, ...
                                            sampleSizeX,...
                                            sampleSizeY,...
                                            numSpinsX,...
                                            numSpinsY,...
                                            T1, ...
                                            T2, ...
                                            chemShift, ...
                                            suppressPlots)
% Calculates and displays the 2D (x & y) frequency reponse of a pulse
% for a given sample, by solving the Bloch equations with relaxation.
% Note that it is assumed the pulses have gradients included, otherwise
% all spatial positions will have the same response!
%
% Inputs:
%
% Var. Name      Units   Description
% pulse          -       Input pulse structure
% zCoordinate    mm      z-position of the spins
% sampleSizeX    mm      Phantom size along x
% sampleSizeY    mm      Phantom size along y
% numSpinsX      -       Num. of spins along x
% numSpinsY      -       Num. of spins along y
% T1             ms      Longitudinal relaxation
% T2             ms      Transverse relaxation
% M0             a.u.    Equilibrium magnetization
% offset         kHz     Chemical shift of spins
% suppressPlots  Boolean Plot results or not?
%
% Outputs: none.

% *********************************
% Create homogeneous square phantom
% *********************************

equilibriumMag = 1; % a.u.

% Create square, homogeneous phantom in x and y. 
spins = initPhantomSquare(sampleSizeX,...
                          sampleSizeY,...
                          zCoordinate, ...
                          numSpinsX,...
                          numSpinsY,...
                          T1, ...
                          T2, ...
                          equilibriumMag, ...
                          chemShift);

% ***********
% Apply pulse
% ***********

spins = ApplyPulse(spins,pulse);


% ****************************
% Extract output magnetization
% ****************************

vecPosX = linspace(-sampleSizeX/2, sampleSizeX/2, numSpinsX);
vecPosY = linspace(-sampleSizeY/2, sampleSizeY/2, numSpinsY);

% Mimic the creation order of the initPhantomSquare routine
% (what I should probably do is use the r-coordinates of the
% spin and interpolate, but I'm too lazy)
curSpin = 0;
for curPosX=1:numSpinsX
    for curPosY=1:numSpinsY
        curSpin = curSpin + 1;
        Mx(curPosX,curPosY) = spins(curSpin).M(1);
        My(curPosX,curPosY) = spins(curSpin).M(2);
        Mz(curPosX,curPosY) = spins(curSpin).M(3);
    end
end

% ************
% Plot results
% ************

if (suppressPlots~=1)
    figure

    subplot(2,2,1)
    imagesc(vecPosX, vecPosY, Mz, [-1 1]);
    colorbar
    title('Mz (a.u.)');
    xlabel('x (mm)');
    ylabel('y (mm)');

    subplot(2,2,2)
    imagesc(vecPosX, vecPosY, Mx, [-1 1]);
    colorbar
    title('Mx (a.u.)');
    xlabel('x (mm)');
    ylabel('y (mm)');

    subplot(2,2,3)
    imagesc(vecPosX, vecPosY, My, [-1 1]);
    colorbar
    title('My (a.u.)');
    xlabel('x (mm)');
    ylabel('y (mm)');

    subplot(2,2,4)
    imagesc(vecPosX, vecPosY, abs(Mx+1i*My), [-1 1]);
    colorbar
    title('|M_{xy}| (a.u.)');
    xlabel('x (mm)');
    ylabel('y (mm)');
end