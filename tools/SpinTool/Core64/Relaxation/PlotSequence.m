function PlotSequence(seq, varargin)
% PlotSequence  Draws given sequence
%   PlotSequence(seq)  Creates a new figure and plots the supplied 
%   sequence.
%
%   PlotSequence(seq, 'propertyname', propertyvalue, ...)  Allows the user 
%   to further supply plotting parameters. These are (case-insensitive):
%     'Figure'     Specifies handle to existing figure. If omitted, a new
%                  figure will be created.
%     'Plot'       A cell array of strings, specifying what to plot. Options
%                  are 'Gx', 'Gy', 'Gz', 'B1Amp', 'B1Freq', 'B1Phase'. For example,
%                  {'B1Amp', 'B1Phase', 'Gx'}. If omitted, only RF amp is 
%                  plotted.
%     'Range'      A 1x2 vector indicating the initial and final plotting 
%                  limits, in temporal units. If omitted, everything is 
%                  plotted.
%     'PulseNames' A cell array of strings of the names of all pulses
%                  in the sequence. An empty string or empty set mean 
%                  the default name will be used for that pulse.
%     'HardPulseDuration'
%   Additional plotting parameter pairs are described in GetSequenceDiagram

p = inputParser;

p.addParameter('HardPulseDuration', 1, @(x) isa(x, 'double'));  
p.addParameter('HardPulseScaling', 1, @(x) isa(x, 'double'));
p.addParameter('HardPulse90Amp', 0, @(x) isa(x, 'double'));
p.addParameter('AMPulsesZeroPhase', true, @(x) isa(x, 'logical'));
p.addParameter('KillMxyDuration', 0, @(x) isa(x, 'double'));
p.addParameter('PurgeRampTime', 0, @(x) isa(x,'double'));
p.addParameter('Accuracy', 0, @(x) isnumeric(x) && (x>=0));
p.addParameter('ShowPulsePhase', true, @(x) islogical(x));
p.addParameter('Figure', []);
p.addParameter('Plot', {'B1Amp'}, @(x) iscell(x));
p.addParameter('Range', [], @(x) isa(x, 'double'));
p.addParameter('PulseNames', []);

p.parse(varargin{:});


if isempty(p.Results.Figure)
    hFig = figure;
else
    hFig = p.Results.Figure;
end
set(hFig, 'Color', [1 1 1]);

[tAxis, B1Amp, B1Phase, Gx, Gy, Gz, pulseCenters, pulseNames, GxSS, GySS, GzSS, GxPurge, GyPurge, GzPurge] = ...
    GetSequenceDiagram(seq, 'HardPulseDuration', p.Results.HardPulseDuration, ...
                            'HardPulseScaling', p.Results.HardPulseScaling, ...
                            'HardPulse90Amp', p.Results.HardPulse90Amp,...
                            'AMPulsesZeroPhase', p.Results.AMPulsesZeroPhase,...
                            'KillMxyDuration', p.Results.KillMxyDuration, ...
                            'PurgeRampTime', p.Results.PurgeRampTime, ...
                            'Accuracy', p.Results.Accuracy, ...
                            'ShowPulsePhase', p.Results.ShowPulsePhase);

if ~isempty(p.Results.PulseNames) && iscell(p.Results.PulseNames)
    for idx=1:numel(p.Results.PulseNames)
        curName = p.Results.PulseNames{idx};
        if ~isempty(curName) 
            pulseNames{idx} = curName;
        end
    end
end
                        
% for idx=1:
% pulseNames = {'90°', '', '180°', '180°', '180°', '180°'};

SSColor = [0 0 0];
SSAlpha = 0.1;
purgeColor = [0 0 0];
purgeAlpha = 0.7;
pulseColor = [ 0 0 0];
pulseAlpha = 0.8;



gradLimits = [min([Gx, Gy, Gz]), max([Gx, Gy, Gz])];
plotVec = p.Results.Plot; 
plotVec = plotVec(numel(plotVec):-1:1);

numPlots = numel(plotVec); 
% plotRelativeHeights = [0.6 0.6 0.6 0.6 1];
plotRelativeHeights = ones(1, numPlots);
plotYSize = plotRelativeHeights/sum(plotRelativeHeights);
plotYHeight = plotYSize*0.8;
plotDeltaY = plotYSize - plotYHeight;
for idx=1:numPlots
    hAxes(idx) = axes('Parent', hFig, 'Position', [0.1 sum(plotYSize(1:idx-1))+plotDeltaY(idx)/2 0.9 plotYHeight(idx)], 'Units', 'normalized', 'YTick', []);
    switch lower(plotVec{idx})
        case 'b1amp'
            x = tAxis(1:end-1);
            y = B1Amp(2:end);
            dY = max(y) - min(y);
            yCenter = (max(y) + min(y))/2;
            yLimits = yCenter + [-dY dY]*1.2/2;
            varName = '^1H ';
            line(x, y, 'Color', [0 0 0], 'LineWidth', 1);
            area(x,y,  'FaceColor', pulseColor, 'FaceAlpha', pulseAlpha);
            % patch(x, y, [0 0 0], 'FaceAlpha', 0.3);
        case 'b1phase'
            x = tAxis(1:end-1);
            y = B1Phase(2:end);
            dY = max(y) - min(y);
            yCenter = (max(y) + min(y))/2;
            yLimits = yCenter + [-dY dY]*1.2/2;
            varName = '\phi_{RF}';
            line(x, y, 'Color', [0 0 0], 'LineWidth', 1);
        case 'b1freq'
            x = tAxis(1:end-1);
            dt = x(2) - x(1);
            y = diff(B1Phase)./dt;
            dY = max(y) - min(y);
            yCenter = (max(y) + min(y))/2;
            yLimits = yCenter + [-dY dY]*1.2/2;
            varName = 'Freq.';
            line(x, y, 'Color', [0 0 0], 'LineWidth', 1);
        case 'gx'
            x = tAxis(1:end-1);
            y = Gx(2:end);
            dY = diff(gradLimits);
            yCenter = mean(gradLimits);
            yLimits = yCenter + [-dY dY]*1.2/2;
            varName = 'G_{RO} ';
            hold(hAxes(idx), 'on');
            line(x, y, 'Color', [0 0 0], 'LineWidth', 1);
            area(x, GxSS(2:end), 'FaceColor', SSColor, 'FaceAlpha', SSAlpha);
            area(x, GxPurge(2:end), 'FaceColor', purgeColor, 'FaceAlpha', purgeAlpha);
            hold(hAxes(idx), 'off');
        case 'gy'
            x = tAxis(1:end-1);
            y = Gy(2:end);
            dY = diff(gradLimits);
            yCenter = mean(gradLimits);
            yLimits = yCenter + [-dY dY]*1.2/2;
            varName = 'G_{PE} ';
            hold(hAxes(idx), 'on');
            line(x, y, 'Color', [0 0 0], 'LineWidth', 1);
            area(x, GySS(2:end), 'FaceColor', SSColor, 'FaceAlpha', SSAlpha);
            area(x, GyPurge(2:end), 'FaceColor', purgeColor, 'FaceAlpha', purgeAlpha);
            hold(hAxes(idx), 'off');
        case 'gz'
            x = tAxis(1:end-1);
            y = Gz(2:end);
            dY = diff(gradLimits);
            yCenter = mean(gradLimits);
            yLimits = yCenter + [-dY dY]*1.2/2;
            varName = 'G_{SL} ';
            hold(hAxes(idx), 'on');
            line(x, y, 'Color', [0 0 0], 'LineWidth', 1);
            area(x, GzSS(2:end), 'FaceColor', SSColor, 'FaceAlpha', SSAlpha);
            area(x, GzPurge(2:end), 'FaceColor', purgeColor, 'FaceAlpha', purgeAlpha);
            hold(hAxes(idx), 'off');
        otherwise
            x = tAxis;
            y = tAxis.*0;
            varName = '';
            line(x, y, 'Color', [0 0 0], 'LineWidth', 1);
    end
    if yLimits(1)<yLimits(2)
        set(hAxes(idx), 'YLim', yLimits);
    end
    if idx>1
        set(hAxes(idx), 'XTick', []);
    else
        
    end
    if strcmpi(lower(plotVec{idx}), 'B1Amp') && ~isempty(pulseCenters)
        for idxLabel=1:numel(pulseNames)
            text(pulseCenters(idxLabel), yLimits(2)*1.05, pulseNames{idxLabel}, 'HorizontalAlignment', 'center');
        end
    end
    set(hAxes(idx), 'XAxisLocation', 'origin', 'YTick', [], 'Box', 'off', 'YColor', [1 1 1], 'XTick', []);
    text(0, 0, varName, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman', 'FontSize', 12);
    
end
linkaxes(hAxes, 'x');
                        