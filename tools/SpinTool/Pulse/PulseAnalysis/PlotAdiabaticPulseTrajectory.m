function [B, M] = PlotAdiabaticPulseTrajectory(pulse, offset, Mi, frame, isPlot, numAnnotatePts)
% Plots the trajectory of the magnetic field in the chosen frame. 
%
% Input Variables
% Variable Name   Description
% pulse           Input pulse object
% offset          Offset of spin (constant), in kHz
% Mi              The magnetization at the beginning of the pulse.
%                 For example, [0; 0; 1] for magnetization along z.
% frame           A (case-insensitive) string:
%                 'rot':    Rotating frame.
%                 'FM':     Frequency modulated frame
%                 'tilted': Tilted frame
%                 'pulse':  Pulse's frame of reference
% isPlot          Optional. 0 or 1 (default). Set to 0 to not plot.
% numAnnotatePts  Number of annotation points (i.e. text boxes). These
%                 points are EQUISPACED in time. Set to 0 (default) to
%                 suppress annotations.

if nargin<6
    numAnnotatePts = 0;
end

if nargin<5
    isPlot = 1;
end

if nargin<4
    frame = 'FM';
end

if nargin<3
    Mi = [0; 0; 1];
end

if nargin<2
    offset = 0;
end

% Calculate the field components, in kHz
numSteps = numel(pulse.RFamp);
dt = pulse.tp/numSteps;
switch lower(frame)
    case 'rot'
        Bx = pulse.RFamp.*cos(pulse.RFphase);
        By = pulse.RFamp.*sin(pulse.RFphase);
        Bz = offset*ones(1,numSteps);
        camAZ = -37.5; % Standard 3D view
        camEL = 30; % Standard 3D view
    case 'fm'
        vRF = diff(pulse.RFphase)/(2*pi*dt);
        vRF(end+1) = vRF(end) + (vRF(end)-vRF(end-1)); % Interpolate last pt
        Bx = pulse.RFamp;
        By = zeros(1,numSteps);
        Bz = offset - vRF;
        camAZ = 0; % View from "side" (i.e. view xz plane)
        camEL = 0;
    case 'tilted'
        vRF = diff(pulse.RFphase)/(2*pi*dt);
        vRF(end+1) = vRF(end) + (vRF(end)-vRF(end-1)); % Interpolate last pt
        BxFM = pulse.RFamp;
        BzFM = offset - vRF;
        tiltAngle = atan2(BxFM, -BzFM);
        dTiltAngle = diff(tiltAngle)/dt;
        dTiltAngle(end+1) = dTiltAngle(end) + (dTiltAngle(end)-dTiltAngle(end-1));
        Bx = dTiltAngle/(2*pi);
        By = zeros(1,numSteps);
        Bz = sqrt(BxFM.^2 + BzFM.^2);
        camAZ = 0; % View from "side" (i.e. view xz plane)
        camEL = 0;
end

% Simulate magnetization's trajectory as a function of time
[Mx, My, Mz] = ApplyFieldDiagnostics(Bx, By, Bz+offset, Mi, pulse.tp);

M = [Mx; My; Mz];
B = [Bx; By; Bz];

% If B is constant, add a point at (0,0,0) so that we can see B
if (sum(abs(diff(Bx)))+ sum(abs(diff(By)))+sum(abs(diff(Bz))))==0 % B is constant
    disp('adding zero');
    Bx = [0, Bx];
    By = [0, By];
    Bz = [0, Bz];
end

BNorm = sqrt(Bx.^2 + By.^2 + Bz.^2);
BNormMax = max(BNorm);
BxNorm = Bx./BNormMax;
ByNorm = By./BNormMax;
BzNorm = Bz./BNormMax;

annotatedPts = round(linspace(1,numSteps,numAnnotatePts));

% Plot resulting trajectory and resulting RF's trajectory
if (isPlot)
    figure
    subplot(1,2,1)
    line(Bx,By,Bz, 'Color', 'r');
    hold on
    line(Bx(annotatedPts), By(annotatedPts), Bz(annotatedPts), 'MarkerEdgeColor', 'none', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    xlabel('Bx (kHz)');
    ylabel('By (kHz)');
    zlabel('Bz (kHz)');
    title('RF Trajectory');
    % set(gca,'DataAspectRatio', [1 1 1]);
    view(camAZ, camEL);

    subplot(1,2,2)
    line(BxNorm, ByNorm, BzNorm, 'Color', 'r');
    hold on
    line(Mx,My,Mz);
    line(Mx(annotatedPts), My(annotatedPts), Mz(annotatedPts), 'MarkerEdgeColor', 'none', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
    line(BxNorm(annotatedPts), ByNorm(annotatedPts), BzNorm(annotatedPts), 'MarkerEdgeColor', 'none', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    title('Magnetization Trajectory (blue)');
    view(camAZ, camEL);
    xlabel('Mx (au)');
    ylabel('My (au)');
    zlabel('Mz (au)');
end