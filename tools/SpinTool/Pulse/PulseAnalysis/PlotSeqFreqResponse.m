function [response, xAxis, hAxes] = PlotSeqFreqResponse(sequences, initMag, plotMin, ...
                                                 plotMax, numPoints, whatToPlot, ...
                                                 grad, chemShiftVec, B1Scaling, ...
                                                 isPlot, seqNames, displayLimits, ...
                                                 scaleFactor, T1, T2, ppmFactor)
% Plots the response of a sequence.
%   PlotSeqFreqResponse(seq)  Plots the frequency response of a sequence.
%   By default this plots Mz, Mx and My as a function of offset from -10
%   to +10 kHz for 501 spins at isocenter. seq can be a sequence or a
%   cell array of sequences, in which case the results will be overlaid
%   on top of each other when possible.
%   
%   PlotSeqFreqResponse(seq, 'name', value) Name-value pairs allow the
%   user to tailor the appearance of the plot, simulation and spin system.
%   See below for a list of potential pairs.
%
%   [response, xAxis, hAxes] = PlotSeqFreqResponse(seq, 'name', value)
%   Returns the simulated responses - i.e., the spins' magnetization
%   vectors for the simulated profiles - in the cell array response.
%
%
%   See also: ApplySequence (on how sequences should be formattedD)


                                             % SYNTAX: [response, xAxis] = PlotSeqFreqResponse
%              (sequences, initMag, plotMin, ...
%               plotMax, numPoints, whatToPlot, ...
%               grad, chemShiftVec, B1Scaling, ...
%               isPlot, seqNames, displayLimits, ...
%               scaleFactor, T1, T2, ppmFactor)
%
% Plots the frequency reponse of one or several sequences, overlaid. This
% assumes infinite T2, T1.
%
% Input
% Variable Name     Units     Description
% sequences         -         A 1xN cell, containing the different seqs.
% initMag           -         Initial magnetiation vector, 3x1
% plotMin, plotMax  kHz, mm   Plot range. If no gradient is specified,
%                             this is taken to be in kHz. Otherwise,
%                             this is taken in mm.
% numPoints         -         Number of points in plot (= num. of offsets)
% whatToPlot        -         'Mz', 'Mx', 'My', 'phase', 'Mxy', 'flipangle'
%                             (not case sensitive)
% grad              'x', 'y'  Optional. If present, the appropriate axis
%                   'z'       gradient will be used, and the 
%                             x-axis will be displayed in mm, ranging
%                             from plotMin to plotMax. Set to [] to
%                             disregard.
% chemShiftVec         kHz    Optional. If a gradient is supplied, this
%                             specifies the constant chemical shift of the
%                             spins. Only meaningful if a gradient is 
%                             present.
% B1Scaling         -         Optional. B1 scaling factor (Representing RF
%                             inhomogeneity). Set to 1 for no
%                             inhomogeneity.
% isPlot            0, 1,     Optional. Set to 0 to suppress plotting,
%                   2, 3      1 to overlay, 2 to juxtapose horizontally,
%                             3 to juxtapose vertically. Default: 1.
% seqNames          1xN cell  Optional. Contains names of sequences.
% displayLimits     kHz, mm   Optional. [] by default. If present,
%                             vertical lines will be inserted into the
%                             plots at x=displayLimits(1) and
%                             x=displayLimits(2)
% scaleFactor       -         Scales vertical limits on graphs (1 by
%                             default; optional)
% T1, T2            ms        Relaxation constants
% ppmFactor         Hz/ppm    If empty or set to 0, nothing will happen.
%                             If a number > 0, it will be treated as a
%                             scaling factor and the freq. axis will be
%                             set to ppm (if no gradient is on). The
%                             ppmFactor will then equal the conversion
%                             factor between Hz<-->ppm. 4.7 ppm will be 
%                             assigned to 0 Hz.

% If it's just a single sequence, put it into a single-cell array

xAxis = linspace(plotMin, plotMax, numPoints); % Display axis: in kHz/ppm (if grad=[]) or mm (if grad='x', 'y' or 'z')

if (~iscell(sequences))
    fprintf('Problem with PlotSeqFreqResponse: input sequence is not a cell array!\n');
    return
end
% If the first element in the sequence is legitimate then we have a single
% sequence, and therefore we're going to put is as a single member in a
% cell array of sequences. 
if ~strcmpi(ReturnSequenceElementType(sequences, 1), 'unknown')
    sequences = {sequences};
end


numSequences = numel(sequences);

if (nargin<7),  grad = []; end
if (isempty(grad))
    % If a gradient is not supplied, the plot limits are assumed to 
    % represent chemical shifts (offsets), EVEN IF a list of chemical
    % shifts are supplied (i.e. if chemShiftVec is supplied)
    chemShiftVec = xAxis; 
else
    if nargin<8
        % No chemical shift is supplied: assume zero chemical shift
        chemShiftVec = 0;
    else
        % Use the given chemical shift vector (i.e. do nothing)
    end
end


if (nargin<9),  B1Scaling = 1; end
if (nargin<10), isPlot = 1; end
if (nargin<11), seqNames = []; end
if (nargin<12), displayLimits = []; end
if (nargin<13), scaleFactor = 1; end
if (nargin<14), T1=1e6; end
if (nargin<15), T2=1e6; end
if (nargin<16), ppmFactor = 0; end

numChemShifts = numel(chemShiftVec);

if ~isempty(displayLimits)
    if size(displayLimits,1) == 1
        displayLimits = repmat(displayLimits, numSequences, 1);
    end
end

if ~iscell(whatToPlot), whatToPlot = {whatToPlot}; end
if ~iscellstr(whatToPlot), error('Unrecognized plot type (should be mz, mx, mxy, etc ... ).'); end

if isempty(seqNames)
    for idxSeq=1:numSequences
        seqNames{idxSeq} = sprintf('Seq%d',idxSeq);
    end
end
if ~isempty(grad)
    switch isPlot
        case {1,2,3}
            for idxCS=1:numChemShifts
                for idxSeq=1:numSequences
                    plotLegend{idxCS, idxSeq} = sprintf('%s (cs=%.2f kHz)', seqNames{idxSeq}, chemShiftVec(idxCS));
                end
            end
    end
else % No gradient
    for idxSeq=1:numSequences
        plotLegend{1, idxSeq} = sprintf('%s', seqNames{idxSeq});
    end
end

eqMag = 1;
if isempty(grad) % Simulate in kHz
    sampleSize = 0;
    sampleOffset = 0;
    sumDim = 1; % Plot the chemical shift dimension in Mxy
    posVec = 0;
else % Gradient is ignored / not supplied. 
    sampleSize = plotMax-plotMin; % mm
    sampleOffset = (plotMax+plotMin)/2; % mm
    sumDim = 2; % Plot the position dimension in Mxy
    posVec = linspace(sampleOffset - sampleSize/2, sampleOffset + sampleSize/2, numPoints);
end
XDir = 'normal';
numPos = numel(posVec);
if isempty(grad)
    spins = InitSpinsRelax(chemShiftVec, 1, sampleSize, initMag, T1, T2, eqMag, sampleOffset, B1Scaling, grad);
else
    spins = InitSpinsRelax(chemShiftVec, numPoints, sampleSize, initMag, T1, T2, eqMag, sampleOffset, B1Scaling, grad);
end

% Initialize for speed
Mx{numSequences}  = zeros(numChemShifts, numPos); 
My{numSequences}  = zeros(numChemShifts, numPos); 
Mz{numSequences}  = zeros(numChemShifts, numPos); 
Mxy{numSequences} = zeros(numChemShifts, numPos); 

% Extract magnetization
for idxSequence=1:numSequences
    spinsOut = ApplySequence(spins, sequences{idxSequence});
    counter = 0;
    for idxCS=1:numChemShifts
        for idxPos=1:numPos
            counter = counter + 1;
            Mx{idxSequence}(idxCS, idxPos)  = spinsOut(counter).M(1);
            My{idxSequence}(idxCS, idxPos)  = spinsOut(counter).M(2);
            Mz{idxSequence}(idxCS, idxPos)  = spinsOut(counter).M(3);
            Mxy{idxSequence}(idxCS, idxPos) = spinsOut(counter).M(1)+1i*spinsOut(counter).M(2);
        end
    end
    % There are numPoints spins in the sample between plotMin and plotMax
    % Hence there are numPoints/(plotMax-plotMin) spins per mm.
    spinsPerMillimeter = numPoints/(plotMax-plotMin);
    cumSumMxy{idxSequence} = abs(cumsum(Mxy{idxSequence}, sumDim))/spinsPerMillimeter;
    maxCumSum(idxSequence) = max(cumSumMxy{idxSequence}(:));
    voxelProfile{idxSequence} = diff(cumSumMxy{idxSequence}, 1, sumDim);
%     voxelProfile
%     voxelProfile{idxSequence}(end+1) = 0;
    maxVoxelProfile(idxSequence) = max(voxelProfile{idxSequence}(:));
    minVoxelProfile(idxSequence) = min(voxelProfile{idxSequence}(:));
    magMxy{idxSequence} = abs(Mxy{idxSequence});
    flipAngle{idxSequence} = acos(Mz{idxSequence})/pi*180;
    phaseMxy{idxSequence} = unwrap(angle(Mxy{idxSequence}), [], sumDim);
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

plotTypes = {'-', '--', '-.', '-.-'}; % For plotting different chemical shifts
if isPlot, figure; end

% Start by assuming no juxtapositioning
numPlotTypes = numel(whatToPlot); 
switch (isPlot)
    case 1 % Overlay
        numHorPlots = 1;
        numVerPlots = numPlotTypes;
        numPlots = numSequences*numPlotTypes;
    case 2 % Juxtapose horizontally
        numHorPlots = numSequences;
        numVerPlots = numPlotTypes;
        numPlots = numVerPlots*numHorPlots;
    case 3 % Juxtapose vertically
        numHorPlots = numPlotTypes;
        numVerPlots = numSequences;
        numPlots = numVerPlots*numHorPlots;
    otherwise
        numPlots = 0;
end

for idxPlotType=1:numPlotTypes  % mz, mx, my, angle, etc ...
    for idxSequence=1:numSequences 
        switch lower(whatToPlot{idxPlotType})
            case 'mz'
                response{idxPlotType, idxSequence} = Mz{idxSequence};
            case 'mx'
                response{idxPlotType, idxSequence} = Mx{idxSequence};
            case 'my'
                response{idxPlotType, idxSequence} = My{idxSequence};
            case 'mxy'
                response{idxPlotType, idxSequence} = magMxy{idxSequence};
            case 'phase'
                response{idxPlotType, idxSequence} = phaseMxy{idxSequence};
            case 'atan'
                response{idxPlotType, idxSequence} = atan2(My{idxSequence}, Mx{idxSequence});
            case 'flipangle'
                response{idxPlotType, idxSequence} = flipAngle{idxSequence};
            case 'summxy'
                response{idxPlotType, idxSequence} = cumSumMxy{idxSequence};
            case 'voxelprofile'
                response{idxPlotType, idxSequence} = voxelProfile{idxSequence};
        end
    end
end

% Example: for two sequences, and 4 plot types:
% idxPlot goes: 1 2 3 4 5 6 7 8
% Pulse number: 1 2 1 2 1 2 1 2
% Plot type #:  1 1 2 2 3 3 4 4
for idxPlot=1:numPlots % = [# vertical plots]*[# horizontal plots]
    idxSequence = mod(idxPlot-1, numSequences)+1;
    idxPlotType = ceil(idxPlot/numSequences);
    switch (isPlot)
        case 1 % Overlay
            curPlotIdx = idxPlotType;
        case 2 % Juxtapose horizontally
            curPlotIdx = idxPlot; 
        case 3 % Juxtapose vertically
            tempIdx = (idxPlot-1)*numHorPlots+1;
            curPlotIdx =  mod(tempIdx, numPlots) + floor(tempIdx/numPlots);
    end
    hAxes(idxPlot) = subplot(numVerPlots, numHorPlots, curPlotIdx);
    set(hAxes(idxPlot), 'XDir', XDir);
    hold on;
    % plot(xAxis, response{idxPlotType, idxSequence},'Color',colorOrder(idxCS, idxSequence,:));
    if isPlot==1 && ~isempty(grad)
        if numChemShifts>1
            plot(xAxis, response{idxPlotType, idxSequence}, 'LineStyle', plotTypes{idxSequence}); % Use different linestyles to differentiate between sequence and chemical shift
        else
            plot(xAxis, response{idxPlotType, idxSequence}, 'Color', colorOrder(idxSequence, :)); % Use different linestyles to differentiate between sequence and chemical shift
        end
    else
        plot(xAxis, response{idxPlotType, idxSequence});
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
            yi = 0; yf = 1*scaleFactor;
            titleStr='|M_{xy}|';
        case 'phase'
             yi = min(phaseMxy{idxSequence})*scaleFactor; yf = max(phaseMxy{idxSequence})*scaleFactor;
            titleStr='Phase of M_{xy}';
        case 'atan'
            titleStr='atan(My/Mx)';
            yi = -pi/2*scaleFactor; yf = pi/2*scaleFactor;
        case 'flipangle'
            titleStr='Flip angle (deg.), modulu pi';
            yi = 0; yf = 180*scaleFactor;
        case 'summxy'
            yi = 0; yf = max(maxCumSum)*1.1*scaleFactor;
            titleStr='|\int M_{xy}|';
        case 'voxelprofile'
            yi = min(minVoxelProfile)*1.1*scaleFactor; yf = max(maxVoxelProfile)*1.1*scaleFactor;
            titleStr='Voxel profile';
    end
    title(titleStr);
    axis([xAxis(1) xAxis(end) yi yf]);
    if (nargin>6) && (~isempty(grad)) 
        xlabel('mm');
    else
        if ppmFactor==0
            xlabel('kHz');
        else
            xlabel('ppm');
        end
    end
    % Add a legend
    switch (isPlot)
        case 1
            legend(plotLegend{:});
        case {2,3}
            legend(plotLegend{:, idxSequence});
    end
    
    if ~isempty(displayLimits)
        plot([displayLimits(idxSequence, 1),displayLimits(idxSequence, 1)], [yi, yf], 'k--');
        plot([displayLimits(idxSequence, 2),displayLimits(idxSequence, 2)], [yi, yf], 'k--');
    end
end

if numPlots>0
    linkaxes(hAxes, 'x');
end
