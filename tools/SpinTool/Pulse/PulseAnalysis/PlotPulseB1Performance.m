function figHandle = PlotPulseB1Performance(pulse, plotType, B1Type, B1Limits, freqLimits, whatToPlot, initMag, figHandle, plotLegend, colorLimits)
% Draws a 2D intensity plot of a given pulse's frequency response for
% different B1 inhomogeneity conditions.
%
% Inputs
%
% Var. Name       Type           Description
% pulse           -              Input pulse.
% plotType        string         '1d', '2d' or '2dcontour'. 
% B1Type          string         'scale' or 'khz'. If 'scale', the vector
%                                B1Limits refers to scaling of the B1
%                                field. If set to 'khz', the vector
%                                B1Limits refers to actual maximal 
%                                amplitude of the B1 field.
% B1Limits        3x1            B1Limits(1) and B1Limits(2) represent
%                                the min. and max. limits B1 (for units,
%                                see B1Type input variable). B1Limits(3) 
%                                is the number of equi-spaced points used.
% freqLimits      1x1 or 3x1     If plotType is set to '1d', freqLimits
%                                is a number for which the magnetization
%                                will be plotted as a function of B1, or
%                                a vector (set) of numbers (frequencies) 
%                                which will then appear as multiple plots.
%                                If plotType is set to '2d', freqLimits(1) 
%                                and freqLimits(2) represent
%                                the min. and max. frequencies of the plot.
%                                freqLimits(3) is the number of equi-spaced 
%                                points used.
% whatToPlot      string         'mz', 'mx', 'my', 'mxy', 'flipangle'
% initMag         3x1            Initial magnetization. Optional. If
%                                omitted, [0; 0; 1] will be assumed.
% figHandle       1x1 real       Optional. Handle to figure. If omitted, a
%                                new figure will be created
% plotLegend      1xnumPulses    Cell array of strings containing the 
%                                legend for the input pulses
% colorLimits     1x2 real       Color limits for figures. If not supplied,
%                                determined automatically. Must be between
%                                0 and 1.

if ispulse(pulse) % A single pulse
    pulse = {pulse};
    numPulses = 1;
    isSeq = 0;
else
    if ispulse(pulse{1}) % A cell array of pulses
        isSeq = 0;
        numPulses = numel(pulse); 
    else % A sequence
        isSeq = 1;
        numPulses = 1; 
    end
end
if numPulses>8
    fprintf('Cannot display more than eight pules. Aborting!');
    return
end
if (nargin<9)
    plotLegend = [];
end
if (nargin<8)
    figHandle = figure;
else
    if isempty(figHandle)
        figHandle = figure;
    end
end
if (nargin<7)
    initMag = [0; 0; 1];
end
if (nargin<6)
    whatToPlot = 'mz';
end    

% Colors in which the different pulses will appear (when plotting 1d).
% This means that no more than 8 pulses can be input
colorOrder = [...
         0         0    1.0000;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500];


for idxPulse=1:numPulses
    if (isSeq==0)
        maxB1(idxPulse) = max(pulse{idxPulse}.RFamp); % kHz
    else
        maxB1(idxPulse) = GetSeqMaxB1(pulse);
    end
    switch lower(B1Type)
        case 'scale'
            B1Vec{idxPulse} = linspace(B1Limits(1), B1Limits(2), B1Limits(3));
        case 'khz'
            % We "cleverly" set the scaling of B1 such that the maximal
            % values of B1 match those specified in B1Limits
            B1VecKhz = linspace(B1Limits(1), B1Limits(2), B1Limits(3));
            B1Vec{idxPulse} = B1VecKhz./maxB1(idxPulse);
        otherwise
            fprintf('Error in PlotPulseB1Performance: input B1Type = %s unrecognized. Aborting! \n', B1Type);
            beep
            return
    end
end

% Calculate the response of the pulse with no B1 inhomogeneity
if ismember(lower(plotType), {'2d', '2dcontour', 'mesh'})
    responseHomogeneous = PlotPulseFreqResponse(pulse, initMag, freqLimits(1), freqLimits(2), 501, whatToPlot, [], 0, 1.0, 0, 0);
    freqAxis1D = linspace(freqLimits(1), freqLimits(2), 501);
end

switch lower(plotType)
    case '1d'
        numFreqs = numel(freqLimits);
        numHorPlots = ceil(sqrt(numFreqs));
        numVerPlots = ceil(numFreqs/numHorPlots);
        figure(figHandle);
        for idxFreq=1:numFreqs
            curFreq = freqLimits(idxFreq);
            subplot(numHorPlots, numVerPlots, idxFreq);
            hold on
            for idxPulse=1:numPulses
                for idxB1=1:numel(B1Vec{idxPulse})
                    spins = InitSpinsRelax(curFreq, 1, 1, initMag, 1e6, 1e6, 1); % Initialize a single spin
                    spins.B1 = B1Vec{idxPulse}(idxB1);
                    if (isSeq==1)
                        [spins, ~] = ApplySequence(spins, pulse);
                    else
                        spins = ApplyPulseRelax(spins, pulse{idxPulse});
                    end
                    switch lower(whatToPlot)
                        case 'mz'
                            M(idxB1) = spins.M(3);
                            ylabel('M_{z}');
                        case 'mx'
                            M(idxB1) = spins.M(1);
                            ylabel('M_{x}');
                        case 'my'
                            M(idxB1) = spins.M(2);
                            ylabel('M_{y}');
                        case 'mxy'
                            M(idxB1) = abs(spins.M(1)+1i*spins.M(2));
                            ylabel('|M_{xy}|');
                    end
                end
                switch (lower(B1Type))
                    case 'scale'
                        plot(B1Vec{idxPulse}, M, 'Color', colorOrder(idxPulse,:));
                        xlabel('B1 Scaling');
                    case 'khz'
                        plot(B1VecKhz, M, 'Color', colorOrder(idxPulse,:));
                        xlabel('kHz');
                end
            end
            title(sprintf('Pulse performance for range of B_1 values, @ %.3f kHz', curFreq));
            if ~isempty(plotLegend)
                legend(plotLegend)
            end
        end
    case {'2d', '2dcontour', 'mesh'}
        for idxPulse=1:numPulses
            freqMin = freqLimits(1);
            freqMax = freqLimits(2);
            numPoints = freqLimits(3);
            freqAxis = linspace(freqMin, freqMax, numPoints);
            magMat{idxPulse} = zeros(numel(B1Vec{idxPulse}), numPoints);
            for idxB1=1:numel(B1Vec{idxPulse})
                spins = InitSpinsRelax(freqAxis, 1, 1, initMag, 1e6, 1e6, 1);
                for idx=1:numPoints
                    spins(idx).B1 = B1Vec{idxPulse}(idxB1);
                end
                if (isSeq==1)
                    [spins, ~] = ApplySequence(spins, pulse);
                else
                    spins = ApplyPulseRelax(spins, pulse{idxPulse});
                end
                for idx=1:numPoints
                    switch lower(whatToPlot)
                        case 'mz'
                            magMat{idxPulse}(idxB1, idx) = spins(idx).M(3);
                            contourVec = [0.6 0.3 0 -0.3 -0.6 -0.9 -0.95 -0.98 -0.99 1];
                            plotTitle = 'Mz';
                        case 'mx'
                            magMat{idxPulse}(idxB1, idx) = spins(idx).M(1);
                            contourVec = 0;
                            plotTitle = 'Mx';
                        case 'my'
                            magMat{idxPulse}(idxB1, idx) = spins(idx).M(2);
                            contourVec = 0;
                            plotTitle = 'My';
                        case 'mxy'
                            magMat{idxPulse}(idxB1, idx) = abs(spins(idx).M(1)+1i*spins(idx).M(2));
                            contourVec = [0 0.01 0.05 0.1 0.4 0.7 0.9 0.95 0.98 0.99 1];
                            plotTitle = '|Mxy|';
                        case 'flipangle'
                            magMat{idxPulse}(idxB1, idx) = acos(spins(idx).M(3))/pi*180;
                            contourVec = [0 1 2 5 10 20 40 70 80 85 88 90 92 95 100 110 120 130 160 170 175 178 179 180];
                            plotTitle = 'Flip angle';
                    end
                end
            end
        end
        globalMax = -1e10; globalMin = 1e10;
        for idxPulse=1:numPulses
            globalMax = max(globalMax, max(magMat{idxPulse}(:)));
            globalMin = min(globalMin, min(magMat{idxPulse}(:)));
        end
end

numRowPlots = 3;
numColPlots = numPulses;

switch lower(whatToPlot)
    case {'mz', 'mx', 'my'}
        histEdges = [-1:0.01:1];
    case 'mxy'
        histEdges = [0:0.01:1];
    case 'flipangle'
        histEdges = [0:1:180];
end    

switch lower(plotType)
    case {'2d', 'mesh'}
        figure(figHandle);
        for idxPulse=1:numPulses
            subplot(numRowPlots, numColPlots, idxPulse);
            plotProfile = magMat{idxPulse};
            response1D = responseHomogeneous{idxPulse};
            switch lower(B1Type)
                case 'scale'
                    switch lower(plotType)
                        case '2d'
                            imagesc(freqAxis, B1Vec{idxPulse}, plotProfile);
                        case 'mesh'
                            mesh(freqAxis, B1Vec{idxPulse}, plotProfile);
                    end
                    ylabel(sprintf('B1 Scaling (Nominal: %.1f kHz)', maxB1(idxPulse)));
                case 'khz'
                    B1VecKhz = linspace(B1Limits(1), B1Limits(2), B1Limits(3));
                    switch lower(plotType)
                        case '2d'
                            imagesc(freqAxis, B1VecKhz, plotProfile);
                        case 'mesh'
                            mesh(freqAxis, B1VecKhz, plotProfile);
                    end
                    ylabel('B1 (kHz)');
            end
            if nargin>=10
                set(gca,'CLim', colorLimits);
            end
            colorbar
            set(gca, 'YDir', 'normal');
            xlabel('Offset (kHz)');
            if ~isempty(plotLegend), title(plotLegend{idxPulse}); end
            if (isSeq==1)
                % title(sprintf('%s. Sequence Duration: %.2f ms. SAR (Rel. to Ref): %.2f.', plotTitle, GetSeqDuration(pulse), GetSeqSAR(pulse, 'ref'))));
                title(plotTitle);
            else
                title(sprintf('%s. Pulse Duration: %.2f ms. SAR (Rel. to Ref): %.2f. Steps: %d. Peak power (No B1 inho): %.2f kHz', plotTitle, pulse{idxPulse}.tp, CalcSAR(pulse{idxPulse}, 'ref'), numel(pulse{idxPulse}.RFamp), max(pulse{idxPulse}.RFamp)));
            end
                
            subplot(numRowPlots, numColPlots, idxPulse + numColPlots);
            [N, BIN] = histc(plotProfile(:), histEdges);
            bar(histEdges, N, 'histc');
            q = quantile(plotProfile(:),[0 .25 .5 .75 1.0 ]);
            title(sprintf('Mean: %.2f, SD: %.2f, Quantiles (0/25/50/75/100): %.2f / %.2f / %.2f / %.2f / %.2f', mean(plotProfile(:)), std(plotProfile(:)), q(1), q(2), q(3), q(4), q(5)));
            xlabel(whatToPlot);
            xlim([globalMin globalMax]);

            subplot(numRowPlots, numColPlots, idxPulse + numColPlots*2);
            plot(freqAxis1D, response1D);
            xlabel('kHz');
            title(sprintf('%s for no B1 inhomogeneity', plotTitle));
            
        end
    case '2dcontour'
        figure(figHandle);
        for idxPulse=1:numPulses
            subplot(numRowPlots, numColPlots, idxPulse);
            plotProfile = magMat{idxPulse};
            response1D = responseHomogeneous{idxPulse};
            switch lower(B1Type)
                case 'scale'
                    if contourVec==0
                        [C, h] = contour(freqAxis, B1Vec{idxPulse}, plotProfile,10);
                    else
                        [C, h] = contour(freqAxis, B1Vec{idxPulse}, plotProfile,contourVec);
                    end
                    clabel(C,h);
                    ylabel({'B1 Scaling', sprintf('(Nominal: %.1f kHz)', maxB1(idxPulse))});
                case 'khz'
                    B1VecKhz = linspace(B1Limits(1), B1Limits(2), B1Limits(3));
                    if contourVec==0
                        [C, h] = contour(freqAxis, B1VecKhz, plotProfile,10);
                    else
                        [C, h] = contour(freqAxis, B1VecKhz, plotProfile,contourVec);
                    end
                    clabel(C,h);
                    ylabel('B1 (kHz)');
            end
            % set(gca, 'CLim', [globalMin globalMax]);
            set(gca, 'YDir', 'normal');
            xlabel('Offset (kHz)');
            if ~isempty(plotLegend), title(plotLegend{idxPulse}); end
            if (isSeq==1)
                % title(sprintf('%s. Sequence Duration: %.2f ms. SAR (Rel. to Ref): %.2f.', plotTitle, GetSeqDuration(pulse), GetSeqSAR(pulse, 'ref'))));
                title(plotTitle);
            else
                title(sprintf('%s. Pulse Duration: %.2f ms. SAR (Rel. to Ref): %.2f. Steps: %d', plotTitle, pulse{idxPulse}.tp, CalcSAR(pulse{idxPulse}, 'ref'), numel(pulse{idxPulse}.RFamp)));
            end
            
            subplot(numRowPlots, numColPlots, idxPulse + numColPlots);
            [N, BIN] = histc(plotProfile(:), histEdges);
            bar(histEdges, N, 'histc');
            q = quantile(plotProfile(:),[0 .25 .5 .75 1.0 ]);
            title(sprintf('Mean: %.2f, SD: %.2f, Quantiles (0/25/50/75/100): %.2f / %.2f / %.2f / %.2f / %.2f', mean(plotProfile(:)), std(plotProfile(:)), q(1), q(2), q(3), q(4), q(5)));
            xlabel(whatToPlot);
            xlim([globalMin globalMax]);
            
            subplot(numRowPlots, numColPlots, idxPulse + numColPlots*2);
            plot(freqAxis1D, response1D);
            xlabel('kHz');
            title(sprintf('%s for no B1 inhomogeneity', plotTitle));
            
        end
end