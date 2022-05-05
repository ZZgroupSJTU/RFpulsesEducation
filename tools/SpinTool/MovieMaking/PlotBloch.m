function [sphereCenter, scaleFactor] = PlotBloch(isInvertedColor, isPlotFID, camPosition, isPlotCircles, hAxes)
% Plots a 3D Bloch Sphere, and returns the handle to the figure.
% If isInvertedColor is set to true, the background will be black and the
% Bloch sphere will be white.
% If isPlotFID is set to true, room will be alotted to a plot of the FID as
% a function of time (as the spins evolve).

if nargin<3, camPosition = []; end


if nargin<4
    isPlotCircles = 1;
end

if isempty(camPosition) 
    camPosition = [0.99 0.1 0.2];
end

if nargin<2
    isPlotFID = 0;
end

if isPlotFID 
    scaleFactor = 0.6;
    sphereCenter = [0 0 0.5];
else
    scaleFactor = 1.0;
    sphereCenter = [0 0 0];
end

if nargin<1
    isInvertedColor = 0;
end

if (isInvertedColor)
    bkColor = [0 0 0];
    fgColor = [1 1 1];
else
    bkColor = [1 1 1];
    fgColor = [0 0 0];
end

if nargin<5, hAxes = gca; end
axes(hAxes);

% Create axes
hold on
if isPlotCircles
    arrowPlot(sphereCenter+[0 0 0], sphereCenter+[0 0 1]*scaleFactor, fgColor, 0.01);
    arrowPlot(sphereCenter+[0 0 0], sphereCenter+[0 1 0]*scaleFactor, fgColor, 0.01);
    arrowPlot(sphereCenter+[0 0 0], sphereCenter+[1 0 0]*scaleFactor, fgColor, 0.01);
end
grid off


% set axes
L = 1.1;
axis([-L L -L L -L L]);

% Set background color
set(hAxes,'Color', bkColor);

if isPlotCircles
    % Create "sphere"
    N=50;
    x=linspace(-1,1,N);
    y1 = sqrt(1-x.^2);
    y2 = -y1;

    x = x.*scaleFactor;
    y1 = y1.*scaleFactor;
    y2 = y2.*scaleFactor;

    plot3(sphereCenter(1)+zeros(1,N), sphereCenter(2)+x, sphereCenter(3)+y1,'Color',fgColor,'LineWidth',1);
    plot3(sphereCenter(1)+zeros(1,N), sphereCenter(2)+x, sphereCenter(3)+y2,'Color',fgColor,'LineWidth',1);

    % Create "great circle"
    N=50;
    x=linspace(-1,1,N);
    y1 = sqrt(1-x.^2);
    y2 = -y1;

    x = x.*scaleFactor;
    y1 = y1.*scaleFactor;
    y2 = y2.*scaleFactor;

    plot3(sphereCenter(1)+x, sphereCenter(2)+y1, sphereCenter(3)+zeros(1,N), 'Color', fgColor, 'LineStyle', ':');
    plot3(sphereCenter(1)+x, sphereCenter(2)+y2, sphereCenter(3)+zeros(1,N), 'Color', fgColor, 'LineStyle', ':');
end


campos(camPosition)