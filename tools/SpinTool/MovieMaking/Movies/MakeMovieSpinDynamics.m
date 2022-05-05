% ========================================================================
% Creates several movies illustrating the dynamics of spins under an
% external magnetic field (no relaxation).
% ========================================================================

clear all
close all
clc

isPlotConstantPrecession = 1;
isPlotMovingPrecession = 1;
isPlotPrettyTrajectory = 0;

% ========================================================================
% Precession around a constant magnetic field
% ========================================================================

if isPlotConstantPrecession
    pulseDuration = 30; % ms
    spinPos = [0;0;0]; % mm
    initialMag = [0.2; 0.2;1];
    numFramesToSkip = 1;
    addAnnotations = 0;
    isNormGlobalB = 0.9;
    isShowNormalPlane = 0;
    isPlotSpinTrajectory = 0;
    isInvertedColors = 1;
    isPlotProgressBar = 0;
    numPulseSteps = 40;
    isPlotBlochSphere = 0;
    B = repmat([0.3; 0.3; 1]*0.031, 1, numPulseSteps);

    outputMovie=MovieBlochB(B, pulseDuration, spinPos,initialMag,numFramesToSkip,addAnnotations, isNormGlobalB, isShowNormalPlane, isPlotSpinTrajectory, isInvertedColors, isPlotBlochSphere, isPlotProgressBar);

    hFig = figure;
    hFig.Color = [0 0 0];
    ha = gca;
    ha.XColor = [0 0 0];
    ha.YColor = [0 0 0];
    ha.ZColor = [0 0 0];
    ha.Box = 'off';
    ha.XTick = [];
    ha.YTick = [];
    ha.ZTick = [];
    ha.YLim = [-0.2 1.5];
    ha.ZLim = [-0.2 1.5];
    filename = 'testmov.gif';
    for N=1:numel(outputMovie)
        image(ha, outputMovie(N).cdata); 
        ha.XColor = [0 0 0];
        ha.YColor = [0 0 0];
        ha.ZColor = [0 0 0];
        ha.Box = 'off';
        ha.XTick = [];
        ha.YTick = [];
        ha.ZTick = [];
        ha.XLim = [200 400];
        ha.YLim = [50 250];
        %ha.ZLim = [200 400];
        colormap(outputMovie(N).colormap); 
        frame = getframe(hFig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if N == 1 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.02); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.02); 
        end         
    end
end
return
%% =======================================================================
% Precession around a moving magnetic field
% ========================================================================

if isPlotMovingPrecession
    pulseDuration = 40; % ms
    spinPos = [0;0;0]; % mm
    initialMag = [0; 0;1];
    numFramesToSkip = 1;
    addAnnotations = 0;
    isNormGlobalB = 0.9;
    isShowNormalPlane = 0;
    isPlotSpinTrajectory = 1;
    magFieldStrength = 1; % 1 = fast, 0.05 = slow
    isInvertedColors = 1;

    numPulseSteps = 250;

    tt = linspace(0,1, numPulseSteps);
    numRevs = 0.5;
    B = [sin(2*pi*tt*numRevs); zeros(1,numPulseSteps); cos(2*pi*tt*numRevs)]*magFieldStrength;

    outputMovie=MovieBlochB(B, pulseDuration, spinPos,initialMag,numFramesToSkip,addAnnotations, isNormGlobalB, isShowNormalPlane, isPlotSpinTrajectory, isInvertedColors);
end

%% =======================================================================
% A pretty trajectory of no scientific merit
% ========================================================================

if isPlotPrettyTrajectory
    pulseDuration = 40; % ms
    spinPos = [0;0;0]; % mm
    initialMag = [0; 0;1];
    initialMag = initialMag./norm(initialMag);
    numFramesToSkip = 1;
    addAnnotations = 0;
    isNormGlobalB = 0.9;
    isShowNormalPlane = 0;
    isPlotSpinTrajectory = 1;
    magFieldStrength = 0.03; % 1 = fast, 0.05 = slow
    isInvertedColors = 1;

    numPulseSteps = 250;

    tt = linspace(0,1, numPulseSteps);
    numRevs = 2.3;
    % B = [zeros(1,numPulseSteps); sin(2*pi*tt*numRevs); ones(1,numPulseSteps)]*magFieldStrength;
    B = [cos(2*pi*tt*numRevs); sin(2*pi*tt*numRevs); sin(2*pi*tt*numRevs)]*magFieldStrength;

    outputMovie=MovieBlochB(B, pulseDuration, spinPos,initialMag,numFramesToSkip,addAnnotations, isNormGlobalB, isShowNormalPlane, isPlotSpinTrajectory, isInvertedColors);
end
