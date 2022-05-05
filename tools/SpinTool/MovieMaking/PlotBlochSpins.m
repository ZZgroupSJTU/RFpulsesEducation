function PlotBlochSpins(isInvertedColor, isPlotFID, camPosition, spins, numPlottedSpins)
% Plots a 3D Bloch Sphere, and returns the handle to the figure.
% If isInvertedColor is set to true, the background will be black and the
% Bloch sphere will be white.
% If isPlotFID is set to true, room will be alotted to a plot of the FID as
% a function of time (as the spins evolve).

if isstruct(spins)
    numSpins = numel(spins);
    M = zeros(3, numSpins);
    for idx=1:numSpins
        M(:,idx) = spins(idx).M';
    end
else
    numSpins = size(spins, 2);
    M = spins;
end

if numPlottedSpins > numSpins
    disp('Too many plotted spins requested. Aborting ... ');
    disp(['(',num2str(numPlottedSpins),' requested, ', num2str(numSpins),' in sample.)']);
    return
end
plottedSpinsVec = floor(linspace(1,numSpins, numPlottedSpins));

M = M./max(sqrt(sum(M(:, plottedSpinsVec).^2,1)));


if isempty(camPosition), camPosition = [0.99 0.1 0.2]; end
if isempty(isPlotFID), isPlotFID = 0; end
if isempty(isInvertedColor), isInvertedColor = 0; end
if isPlotFID, 
    scaleFactor = 0.7;
else
    scaleFactor = 1.0;
end

hh=figure;
PlotBloch(isInvertedColor, isPlotFID, camPosition);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
if isInvertedColor
    set(gcf,'color','black');
    set(gca,'color','black');
else
    set(gcf,'color','white');
    set(gca,'color','white');
end
set(gca,'LineStyleOrder','--');
set(gca,'Visible','off');
QQ = get(hh,'Position');

% Assign different colors to the different plotted spins.
spinColor = [0 0 1];   % Blue
for curPlottedSpin=1:numPlottedSpins
    if isInvertedColor
        spinColorVec{curPlottedSpin} = spinColor +  (1-spinColor)*(numPlottedSpins - curPlottedSpin)/numPlottedSpins;
    else
        spinColorVec{curPlottedSpin} = spinColor - spinColor*(numPlottedSpins - curPlottedSpin)/numPlottedSpins;
    end
end

% Plot spins
for curPlottedSpin=1:numPlottedSpins
    curM = M(:,plottedSpinsVec(curPlottedSpin));
    q1(curPlottedSpin)=arrowPlot([0 0 0],curM',spinColorVec{curPlottedSpin},0.015);
end;

