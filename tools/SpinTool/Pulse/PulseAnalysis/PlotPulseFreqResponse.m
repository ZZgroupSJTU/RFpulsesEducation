function [response, xAxis] = PlotPulseFreqResponse(pulses, initMag, plotMin, ... 
                                          plotMax, numPoints, whatToPlot, ...
                                          grad, chemShift, B1Scaling, ...
                                          phaseAdjust, isPlot, plotLegend, ...
                                          displayLimits, scaleFactor, ...
                                          T1, T2, plotScale)
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

numPulses = numel(pulses);

if (nargin<7)
    grad = [];
end

if (nargin<8)
    chemShift = 0;
end

if (nargin<9)
    B1Scaling = 1;
end

if (nargin<10)
    phaseAdjust = 0;
end

if (nargin<11)
    isPlot = 1;
end

if (nargin<12)
    plotLegend = [];
end

if (nargin<13)
    displayLimits = [];
end

if (nargin<14)
    scaleFactor = 1;
end

if (nargin<15)
    T1=1e6;
end

if (nargin<16)
    T2=1e6;
end

if (nargin<17)
    plotScale='linear';
end

if ~isempty(displayLimits)
    if size(displayLimits,1) == 1
        displayLimits = repmat(displayLimits, numPulses, 1);
    end
end

if ~iscell(whatToPlot)
    whatToPlot = {whatToPlot};
end


eqMag = 1;
if ((nargin<7) || isempty(grad)) % Gradient is not used. All spins at z=x=y=0.
    for idx=1:numPulses
        pulses{idx}.Gz = zeros(1,numel(pulses{idx}.Gz));
        pulses{idx}.Gy = zeros(1,numel(pulses{idx}.Gy));
        pulses{idx}.Gx = zeros(1,numel(pulses{idx}.Gx));
    end
    chemShiftVec = linspace(plotMin, plotMax, numPoints);
    numSpinsPerShift = 1;
    xAxis = linspace(plotMin, plotMax, numPoints); % in kHz
    sampleSize = 0;
else % Gradient is specified
    chemShiftVec = chemShift;
    numSpinsPerShift = numPoints;
    sampleSize = 0;
    xAxis = linspace(plotMin, plotMax, numPoints); % in mm
end
spins = InitSpinsRelax(chemShiftVec, numSpinsPerShift, sampleSize, initMag, T1, T2, eqMag);
for idx=1:numPoints
    spins(idx).B1 = B1Scaling;
end
if ((nargin>=7) && (~isempty(grad))) % Gradient is specified: put equally-spaced spins along spatial axis of choice.
    posVec = linspace(plotMin, plotMax, numPoints);
    switch lower(grad)
        case 'x'
            for idx=1:numPoints
                spins(idx).r = [posVec(idx); 0; 0];
            end
        case 'y'
            for idx=1:numPoints
                spins(idx).r = [0; posVec(idx); 0];
            end
        case 'z'
            for idx=1:numPoints
                spins(idx).r = [0; 0; posVec(idx)];
            end
    end
end

for idxPulse=1:numPulses
    switch phaseAdjust
        case 2 % Crushers
            spinsOut = PurgeMoment(spins, 0, 0, 1500, 0.01);
            spinsOut = ApplyPulseRelax(spinsOut, pulses{idxPulse});
            spinsOut = PurgeMoment(spinsOut, 0, 0, 1500, 0.01);
        otherwise
            spinsOut = ApplyPulseRelax(spins, pulses{idxPulse});
    end
    for idx=1:numPoints
        Mx{idxPulse}(idx)  = spinsOut(idx).M(1);
        My{idxPulse}(idx)  = spinsOut(idx).M(2);
        Mz{idxPulse}(idx)  = spinsOut(idx).M(3);
        Mxy{idxPulse}(idx) = spinsOut(idx).M(1)+1i*spinsOut(idx).M(2);
    end
    if (phaseAdjust==1) && isempty(grad) % Refocusing (wT/2): in frequency domain
        Mxy{idxPulse} = Mxy{idxPulse}.*exp(2*pi*1i*chemShiftVec*pulses{idxPulse}.tp/2);
        Mx{idxPulse} = real(Mxy{idxPulse});
        My{idxPulse} = imag(Mxy{idxPulse});
    end
    if (phaseAdjust==1) && ~isempty(grad) % Refocusing (wT/2): in frequency domain
        switch lower(grad)
            case 'x'
                G = pulses{idxPulse}.Gx(1); % kHz/mm
            case 'y'
                G = pulses{idxPulse}.Gy(1); % kHz/mm
            case 'z'
                G = pulses{idxPulse}.Gz(1); % kHz/mm
        end
        Mxy{idxPulse} = Mxy{idxPulse}.*exp(2*pi*1i*posVec*G*pulses{idxPulse}.tp/2);
        Mx{idxPulse} = real(Mxy{idxPulse});
        My{idxPulse} = imag(Mxy{idxPulse});
    end
    if (phaseAdjust==3) && isempty(grad) % Refocusing (wT/2): in frequency domain
        Mxy{idxPulse} = Mxy{idxPulse}.*exp(2*pi*1i*chemShiftVec*pulses{idxPulse}.tp);
        Mx{idxPulse} = real(Mxy{idxPulse});
        My{idxPulse} = imag(Mxy{idxPulse});
    end
    if (phaseAdjust==3) && ~isempty(grad) % Refocusing (wT/2): in frequency domain
        switch lower(grad)
            case 'x'
                G = pulses{idxPulse}.Gx(1); % kHz/mm
            case 'y'
                G = pulses{idxPulse}.Gy(1); % kHz/mm
            case 'z'
                G = pulses{idxPulse}.Gz(1); % kHz/mm
        end
        Mxy{idxPulse} = Mxy{idxPulse}.*exp(2*pi*1i*posVec*G*pulses{idxPulse}.tp);
        Mx{idxPulse} = real(Mxy{idxPulse});
        My{idxPulse} = imag(Mxy{idxPulse});
    end
    % There are numPoints spins in the sample between plotMin and plotMax
    % Hence there are numPoints/(plotMax-plotMin) spins per mm.
    spinsPerMillimeter = numPoints/(plotMax-plotMin);
    cumSumMxy{idxPulse} = abs(cumsum(Mxy{idxPulse}))/spinsPerMillimeter;
    maxCumSum(idxPulse) = max(cumSumMxy{idxPulse});
    voxelProfile{idxPulse} = diff(cumSumMxy{idxPulse});
    voxelProfile{idxPulse}(end+1) = 0;
    maxVoxelProfile(idxPulse) = max(voxelProfile{idxPulse});
    minVoxelProfile(idxPulse) = min(voxelProfile{idxPulse});
    magMxy{idxPulse} = abs(Mxy{idxPulse});
    flipAngle{idxPulse} = acos(Mz{idxPulse})/pi*180;
    phaseMxy{idxPulse} = phase(Mxy{idxPulse});
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

if isPlot, hFig=figure; end

% Start by assuming no juxtapositioning
numPlotTypes = numel(whatToPlot); 
switch (isPlot)
    case 1 % Overlay
        numHorPlots = 1;
        numVerPlots = numPlotTypes;
        numPlots = numPulses*numPlotTypes;
    case 2 % Juxtapose horizontally
        numHorPlots = numPulses;
        numVerPlots = numPlotTypes;
        numPlots = numVerPlots*numHorPlots;
    case 3 % Juxtapose vertically
        numHorPlots = numPlotTypes;
        numVerPlots = numPulses;
        numPlots = numVerPlots*numHorPlots;
    otherwise
        numPlots = 0;
end

for idxPlotType=1:numPlotTypes
    for idxPulse=1:numPulses
        switch lower(whatToPlot{idxPlotType})
            case 'mz'
                response{idxPlotType, idxPulse} = Mz{idxPulse};
            case 'mx'
                response{idxPlotType, idxPulse} = Mx{idxPulse};
            case 'my'
                response{idxPlotType, idxPulse} = My{idxPulse};
            case 'mxy'
                response{idxPlotType, idxPulse} = magMxy{idxPulse};
            case 'phase'
                response{idxPlotType, idxPulse} = phaseMxy{idxPulse};
            case 'atan'
                response{idxPlotType, idxPulse} = atan2(My{idxPulse}, Mx{idxPulse});
            case 'atandeg'
                response{idxPlotType, idxPulse} = atan2(My{idxPulse}, Mx{idxPulse})/pi*180;
            case 'flipangle'
                response{idxPlotType, idxPulse} = flipAngle{idxPulse};
            case 'summxy'
                response{idxPlotType, idxPulse} = cumSumMxy{idxPulse};
            case 'voxelprofile'
                response{idxPlotType, idxPulse} = voxelProfile{idxPulse};
        end
    end
end


% idxPlot goes: 1 2 3 4 5 6 7 8
% Pulse number: 1 2 1 2 1 2 1 2
% Plot type #:  1 1 2 2 3 3 4 4
if isPlot
    for idxPlot=1:numPlots % = [# vertical plots]*[# horizontal plots]
        idxPulse = mod(idxPlot-1, numPulses)+1;
        idxPlotType = ceil(idxPlot/numPulses);
        switch (isPlot)
            case 1 % Overlay
                curPlotIdx = idxPlotType;
            case 2 % Juxtapose horizontally
                curPlotIdx = idxPlot; 
            case 3 % Juxtapose vertically
                tempIdx = (idxPlot-1)*numHorPlots+1;
                curPlotIdx =  mod(tempIdx, numPlots) + floor(tempIdx/numPlots);
        end
        ax(idxPlot) = subplot(numVerPlots, numHorPlots, curPlotIdx);
        hold on;
        if strcmpi(plotScale, 'log') && strcmpi(whatToPlot{idxPlotType}, 'mxy')
            plot(xAxis, log10(response{idxPlotType, idxPulse}),'Color',colorOrder(idxPulse,:), 'LineWidth', 1.5);
        else
            plot(xAxis, response{idxPlotType, idxPulse},'Color',colorOrder(idxPulse,:), 'LineWidth', 1.5);
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
                if strcmpi(plotScale, 'linear')
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
        axis([xAxis(1) xAxis(end) yi yf]);
        if (nargin>6) && (~isempty(grad)) 
            xlabel('mm');
        else
            xlabel('kHz');
        end
        if ~isempty(plotLegend)
            switch (isPlot)
                case 1
                    legend(plotLegend); 
                case 2
                    legend(plotLegend{idxPulse});
                case 3
                    legend(plotLegend{idxPulse});
            end
        end
        if ~isempty(displayLimits)
            plot([displayLimits(idxPulse, 1),displayLimits(idxPulse, 1)], [yi, yf], 'k--');
            plot([displayLimits(idxPulse, 2),displayLimits(idxPulse, 2)], [yi, yf], 'k--');
        end
    end

    if numPlots>0
        linkaxes(ax, 'x');
    end

    if numPulses==1
        set(gcf, 'Name', sprintf('Duration: %.2f (ms), Peak B1: %.2f (kHz), Num Steps: %d\n', pulses{1}.tp, max(pulses{1}.RFamp), numel(pulses{1}.RFamp)));
    end
    
    hFig.Color = [ 1 1 1];
end

