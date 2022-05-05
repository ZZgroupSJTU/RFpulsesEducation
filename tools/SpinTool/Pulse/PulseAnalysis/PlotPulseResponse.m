function response = PlotPulseResponse(pulses, varargin)

% SYNTAX: function response = PlotPulseFreqResponse(pulses, initMag, plotMin, ... 
%                                           plotMax, numPoints, whatToPlot, ...
%                                           grad, chemShift, B1Scaling, ...
%                                           phaseAdjust, isPlot, plotLegend, ...
%                                           displayLimits, scaleFactor, ...
%                                           T1, T2)
% Plots the frequency reponse of one or several pulses, overlaid. This
% assumes infinite T2, T1.
%
% Input
% Variable Name     Units     Description
% pulses            -         A 1xN cell, containing the different pulses
% initMag           -         Initial magnetiation vector, 3x1
% plotMin, plotMax  kHz, mm   Plot range. If no gradient is specified,
%                             this is taken to be in kHz. Otherwise,
%                             this is taken in mm.
% numPoints         -         Number of points in plot (= num. of offsets)
% whatToPlot        -         'Mz', 'Mx', 'My', 'phase', 'Mxy', 'flipangle'
%                             'atan', 'atandeg', 'summxy', 'voxelprofile'
%                             (case insensitive)
% grad              'x', 'y'  Optional. If present, the appropriate axis
%                   or 'z'    gradient will be used, and the 
%                             x-axis will be displayed in mm, ranging
%                             from plotMin to plotMax.
% chemShift         kHz       Optional. If a gradient is supplied, this
%                             specifies the constant chemical shift.
% B1Scaling         -         Optional. B1 scaling factor (Representing RF
%                             inhomogeneity). Set to 1 for no
%                             inhomogeneity.
% phaseAdjust       -         Optional. 0 by default.
%                             0: do nothing.
%                             1: A linear phase equal to wT/2 will be added
%                                to the transverse magnetization to "refocus"
%                                the phase and allow for better visualization
%                             2: Use crusher gradients before and after
%                                pulse (usually for 180 pulses)
% isPlot            0, 1,     Optional. Set to 0 to suppress plotting,
%                   2, 3      1 to overlay, 2 to juxtapose horizontally,
%                             3 to juxtapose vertically. Default: 1.
% plotLegend        1xN cell  Optional. Contains names of pulses.
% displayLimits     kHz, mm   Optional. [] by default. If present,
%                             vertical lines will be inserted into the
%                             plots at x=displayLimits(1) and
%                             x=displayLimits(2)
% scaleFactor       -         Scales vertical limits on graphs (1 by
%                             default; optional)
% T1, T2            ms        Relaxation constants (if not provided,
%                             default to 1e6)
% plotScaling       'linear'  Whether to show results in linear or 
%                   'log'     log-10 scale. Only makes sense for
%                             positive quantities, i.e. |Mxy|.

% If it's just a single pulse, put it into a single-cell array
if (~iscell(pulses))
    pulses = {pulses};
end

p = inputParser;
p.addParameter('NumPoints', 101, @(x) isnumeric(x));
p.addParameter('Offset', 0, @(x) isvector(x));
p.addParameter('Position', 0, @(x) isvector(x));
p.addParameter('Legend', [], @(x) iscellstr(x));
p.addParameter('T1', 1e6, @(x) isvector(x));
p.addParameter('T2', 1e6, @(x) isvector(x));
p.addParameter('Scale', 'linear', @(x) ismember(lower(x), {'linear', 'log'}));
p.addParameter('Plotting', 'overlay', @(x) ismember(lower(x), {'overlay', 'horizontal', 'vertical', 'none'}));
p.addParameter('GradAxis', [], @(x) ismember(lower(x), {'x', 'y', 'z', []}));
p.addParameter('InitMag', [0;0;1], @(x) isvector(x));
p.addParameter('PlotTarget', 'mz', @(x) iscellstr(x) || ismember(lower(x), {'mz', 'my', 'mx', 'phase', 'mxy', 'flipangle', 'atan', 'atandeg', 'summxy'}));
p.addParameter('PhaseCycle', 'none', @(x) ismember(lower(x), {'none', '0,180'}));
p.addParameter('B1', 1.0 , @(x) isvector(x));
p.addParameter('Axes', [], @(x) ishandle(x));
p.addParameter('ScaleFactor', 1, @(x) isnumeric(x));
p.parse(varargin{:});


numPulses = numel(pulses);

if iscell(p.Results.PlotTarget)
    whatToPlot = p.Results.PlotTarget;
else
    whatToPlot = {p.Results.PlotTarget};
end

initMag = p.Results.InitMag;
scaleFactor = p.Results.ScaleFactor;
if isempty(p.Results.GradAxis), gradAxis = 'z'; else, gradAxis = lower(p.Results.GradAxis); end
offsetVec = p.Results.Offset;
B1Vec = p.Results.B1;
T1Vec = p.Results.T1;
T2Vec = p.Results.T2;
if isempty(p.Results.Position), posVec = 0; else, posVec = p.Results.Position; end
numVectors = (numel(offsetVec)>1) + (numel(B1Vec)>1) + (numel(T1Vec)>1) + (numel(T2Vec)>1) + (numel(posVec)>1);
if numVectors>2
    error('Cannot plot more than 2 different quantities simultaneously.');
end

M0 = 0;
counter = 0;
for idxPos=1:numel(posVec)
    for idxB1=1:numel(B1Vec)
        for idxT1=1:numel(T1Vec)
            for idxT2=1:numel(T2Vec)
                for idxOffset=1:numel(offsetVec)
                    counter = counter + 1;
                    switch gradAxis
                        case 'x'
                            spins(counter).r = [posVec(idxPos); 0; 0];
                        case 'y'
                            spins(counter).r = [0; posVec(idxPos); 0];
                        case 'z'
                            spins(counter).r = [0; 0; posVec(idxPos)];
                    end
                    spins(counter).M  = initMag;
                    spins(counter).cs = offsetVec(idxOffset); 
                    spins(counter).T1 = T1Vec(idxT1);
                    spins(counter).T2 = T2Vec(idxT2);
                    spins(counter).M0 = M0;
                    spins(counter).B1 = B1Vec(idxB1);
                    spins(counter).B0 = 0;
                    spins(counter).RS = 1;
                end
            end
        end
    end
end


for idxPulse=1:numPulses
    switch lower(p.Results.PhaseCycle)
        case 'none'
            seq = {pulses{idxPulse}};
        case '0,180'
            seq = {{'pulse', pulses{idxPulse}, [0 180], [1 1]/2}};
        otherwise
            error('Unrecognized phase cycle %s', p.Results.PhaseCycle);
    end
    spinsOut = ApplySequence(spins, seq);
    counter = 0;
    for idxPos=1:numel(posVec)
        for idxB1=1:numel(B1Vec)
            for idxT1=1:numel(T1Vec)
                for idxT2=1:numel(T2Vec)
                    for idxOffset=1:numel(offsetVec)
                        counter = counter + 1;
                        Mz(idxPos, idxB1, idxT1, idxT2, idxOffset) = spins(counter).M(3);
                        Mxy(idxPos, idxB1, idxT1, idxT2, idxOffset) = spins(counter).M(1) + 1i*spins(counter).M(2);
                    end
                end
            end
        end
    end
    
    % There are numPoints spins in the sample between plotMin and plotMax
    % Hence there are numPoints/(plotMax-plotMin) spins per mm.
    for idxTarget=1:numel(whatToPlot)
        switch lower(whatToPlot{idxTarget})
            case 'mz'
                response{idxTarget,idxPulse} = squeeze(Mz);
            case 'mx'
                response{idxTarget,idxPulse} = squeeze(real(Mxy));
            case 'my'
                response{idxTarget,idxPulse} = squeeze(imag(Mxy));
            case 'mxy'
                response{idxTarget,idxPulse} = squeeze(abs(Mxy));
            case 'flipangle'
                response{idxTarget,idxPulse} = acos(squeeze(Mz))/pi*180;
            otherwise
                error('Unrecognized target %s\n', whatToPlot{idxTarget});
        end
    end
end

% Support up to 7 plots
colorOrder = [...
         0         0    1.0000;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500];


% Start by assuming no juxtapositioning
numPlotTypes = numel(whatToPlot); 
switch lower(p.Results.Plotting)
    case 'overlay'
        numHorPlots = 1;
        numVerPlots = numPlotTypes;
        numPlots = numPulses*numPlotTypes;
    case 'horizontal' % Juxtapose horizontally
        numHorPlots = numPulses;
        numVerPlots = numPlotTypes;
        numPlots = numVerPlots*numHorPlots;
    case 'vertical' % Juxtapose vertically
        numHorPlots = numPlotTypes;
        numVerPlots = numPulses;
        numPlots = numVerPlots*numHorPlots;
    otherwise
        numPlots = 0;
end

if numPlots==0, return; end

figure;

% idxPlot goes: 1 2 3 4 5 6 7 8
% Pulse number: 1 2 1 2 1 2 1 2
% Plot type #:  1 1 2 2 3 3 4 4
for idxPlot=1:numPlots % = [# vertical plots]*[# horizontal plots]
    idxPulse = mod(idxPlot-1, numPulses)+1;
    idxPlotType = ceil(idxPlot/numPulses);
    switch lower(p.Results.Plotting)
        case 'overlay'
            curPlotIdx = idxPlotType;
        case 'horizontal' % Juxtapose horizontally
            curPlotIdx = idxPlot; 
        case 'vertical' % Juxtapose vertically
            tempIdx = (idxPlot-1)*numHorPlots+1;
            curPlotIdx =  mod(tempIdx, numPlots) + floor(tempIdx/numPlots);
    end
    ax(idxPlot) = subplot(numVerPlots, numHorPlots, curPlotIdx);
    hold on;
    if strcmpi(p.Results.Scale, 'log') && strcmpi(whatToPlot{idxPlotType}, 'mxy')
        plot(posVec, log10(response{idxPlotType, idxPulse}),'Color',colorOrder(idxPulse,:));
    else
        plot(posVec, response{idxPlotType, idxPulse},'Color',colorOrder(idxPulse,:));
    end
    switch lower(whatToPlot{idxPlotType})
        case 'mz'
            yi = -1*scaleFactor; yf = 1*scaleFactor;
            titleStr='M_z';
        case 'mx'
            yi = -1*scaleFactor; yf = 1*scaleFactor;
            titleStr='M_x';
        case 'my'
            yi = -1*scaleFactor; yf = 1*scaleFactor;
            titleStr='M_y';
        case 'mxy'
            if strcmpi(p.Results.Scale, 'linear')
                yi = 0; yf = 1*scaleFactor;
            else
                yi = -6; yf = 0;
            end
            titleStr='|M_{xy}|';
        case 'phase'
             yi = min(phaseMxy{idxPulse})*scaleFactor; yf = max(phaseMxy{idxPulse})*scaleFactor;
            titleStr='Phase of M_{xy}';
        case 'atan'
            titleStr='atan(My/Mx)';
            yi = -pi*scaleFactor; yf = pi*scaleFactor;
        case 'atandeg'
            titleStr='atan(My/Mx) (in deg.)';
            yi = -180*scaleFactor; yf = 180*scaleFactor;
        case 'flipangle'
            titleStr='Flip angle (deg.), modulu pi';
            yi = 0; yf = 180*scaleFactor;
        case 'summxy'
            yi = 0; yf = max(maxCumSum)*1.1*scaleFactor;
            titleStr='Integral of abs(Mxy)';
        case 'voxelprofile'
            yi = min(minVoxelProfile)*1.1*scaleFactor; yf = max(maxVoxelProfile)*1.1*scaleFactor;
            titleStr='Voxel profile';
    end
    title(titleStr);
    axis([posVec(1) posVec(end) yi yf]);
    if (nargin>6) && (~isempty(grad)) 
        xlabel('mm');
    else
        xlabel('kHz');
    end
    if ~isempty(p.Results.Legend)
        switch lower(p.Results.Plotting)
            case 'overlay'
                legend(p.Results.Legend); 
            case {'horizontal', 'vertical'}
                legend(p.Results.Legend{idxPulse});
        end
    end
end

if numPlots>0
    linkaxes(ax, 'x');
end

if numPulses==1
    set(gcf, 'Name', sprintf('Duration: %.2f (ms), Peak B1: %.2f (kHz), Num Steps: %d\n', pulses{1}.tp, max(pulses{1}.RFamp), numel(pulses{1}.RFamp)));
end
