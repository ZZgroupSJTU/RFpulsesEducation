function [response, posVec, csVec] = PlotPulseSpatialSpectralResponse(pulses, initMag, freqLimits, spatialLimits, ... 
                                                     whatToPlot, grad, plotType, B1Scaling, ...
                                                     phaseAdjust, T1, T2)
% SYNTAX: 
%
%     response = PlotPulseSpatialSpectralResponse(pulses, initMag, freqLimits, spatialLimits, ... 
%                                                 whatToPlot, grad, B1Scaling, ...
%                                                 phaseAdjust, T1, T2, plotScale)
%
% Plots the frequency reponse of one or several pulses, overlaid. This
% assumes infinite T2, T1.
%
% Input
% Variable Name     Units     Description
% pulses            -         A 1xN cell, containing the different pulses
% initMag           -         Initial magnetiation vector, 3x1
% freqLimits        kHz       A 1x3 vector, such that freqLimits(1) is
%                             the minimal frequency, freqLimits(2) the max
%                             frequency and freqLimits(3) the number of pts
% spatialLimits     mm        A 1x3 vector, such that spatialLimits(1) is
%                             the minimal pos., spatialLimits(2) the max
%                             pos. and spatialLimits(3) the number of pts
% whatToPlot        -         'Mz', 'Mx', 'My', 'phase', 'Mxy', 'flipangle'
%                             'atan', 'atandeg'
%                             (case insensitive)
% grad              'x', 'y'  If set to 'x', 'y', 'z', the appropriate axis
%                   or 'z',   gradient will be used. If a number is 
%                   OR a      supplied, the existing gradient will be 
%                   number    disregarded, and the supplied gradient in 
%                   in mT/m   mT/m will be used instead. If grad is a 
%                             vector, the i-th element will be used for
%                             the i-th pulse.
% plotType          -         Optional. Either '2d', 'mesh' or 'contour' 
%                             (2d by default). Affects appearance of final
%                             image, but nothing else.
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
% T1, T2            ms        Relaxation constants (if not provided,
%                             default to 1e6)

if (nargin<7), plotType = '2d'; end
if (nargin<8), B1Scaling = 1; end
if (nargin<9), phaseAdjust = 0; end
if (nargin<10), T1=1e8; end
if (nargin<11), T2=1e8; end

% If it's just a single pulse, put it into a single-cell array
if (~iscell(pulses)), pulses = {pulses}; end
numPulses = numel(pulses);

if ~iscell(whatToPlot), whatToPlot = {whatToPlot}; end

eqMag = 1;
csVec = linspace(freqLimits(1), freqLimits(2), freqLimits(3));
posVec = linspace(spatialLimits(1), spatialLimits(2), spatialLimits(3));

if isnumeric(grad)
    if numel(grad)~=numPulses
        grad = ones(1,numPulses)*grad;
    end
    for idx=1:numPulses
        N = numel(pulses{idx}.Gx);
        pulses{idx}.Gx = zeros(1, N);
        pulses{idx}.Gy = zeros(1, N);
        pulses{idx}.Gz = ones(1, N).*grad(idx)*GetGyromagneticRatio('1h')/1000;
    end
end
    


counter = 0;
for idxFreq=1:freqLimits(3)
    for idxPos=1:spatialLimits(3)
        counter = counter + 1;
        if ischar(grad)
            switch lower(grad)
                case 'x'
                    spins(counter).r = [posVec(idxPos); 0; 0];
                case 'y'
                    spins(counter).r = [0; posVec(idxPos); 0];
                case 'z'
                    spins(counter).r = [0; 0; posVec(idxPos)];
                otherwise
            end
        else
            spins(counter).r = [0; 0; posVec(idxPos)];
        end
        spins(counter).M  = initMag;
        spins(counter).cs = csVec(idxFreq);  % in kHz!
        spins(counter).T1 = T1; % ms
        spins(counter).T2 = T2; % ms
        spins(counter).M0 = eqMag; % a.u.
        spins(counter).B1 = B1Scaling; % Scales RF
        spins(counter).B0 = 0; % Offset, in kHz
        spins(counter).RS = 1; % Receiver sensitivity
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
    counter = 0;
    for idxFreq=1:freqLimits(3)
        for idxPos=1:spatialLimits(3)
            counter = counter + 1;
            Mx{idxPulse}(idxFreq, idxPos)  = spinsOut(counter).M(1);
            My{idxPulse}(idxFreq, idxPos)  = spinsOut(counter).M(2);
            Mz{idxPulse}(idxFreq, idxPos)  = spinsOut(counter).M(3);
            Mxy{idxPulse}(idxFreq, idxPos) = spinsOut(counter).M(1) + 1i*spinsOut(counter).M(2);
        end
    end
    if (phaseAdjust==1) % Refocusing (wT/2): in frequency domain
        Mxy{idxPulse} = Mxy{idxPulse}.*exp(2*pi*1i*chemShiftVec*pulses{idxPulse}.tp/2);
        Mx{idxPulse} = real(Mxy{idxPulse});
        My{idxPulse} = imag(Mxy{idxPulse});
    end
end

figure
counter = 0;
for idxPulse=1:numPulses
    for idxPlotType=1:numel(whatToPlot)
        counter = counter + 1;
        ax(counter)=subplot(numPulses, numel(whatToPlot), counter);
        switch lower(whatToPlot{idxPlotType})
            case 'mz'
                switch lower(plotType)
                    case '2d'
                        imagesc(posVec, csVec, Mz{idxPulse});
                    case 'mesh'
                        mesh(posVec, csVec, Mz{idxPulse});
                    case 'contour'
                        contour(posVec, csVec, Mz{idxPulse});
                end
                set(ax(counter), 'CLim', [-1 1]);
                title('M_z');
                response{idxPulse, idxPlotType} = Mz{idxPulse};
            case 'flipangle'
                flipangle = acos(Mz{idxPulse})/pi*180;
                switch lower(plotType)
                    case '2d'
                        imagesc(posVec, csVec, flipangle);
                    case 'mesh'
                        mesh(posVec, csVec, flipangle);
                    case 'contour'
                        contour(posVec, csVec, flipangle);
                end
                set(ax(counter), 'CLim', [-1 1]);
                title('M_z');
                response{idxPulse, idxPlotType} = flipangle;
            case 'mx'
                switch lower(plotType)
                    case '2d'
                        imagesc(posVec, csVec, Mx{idxPulse});
                    case 'mesh'
                        mesh(posVec, csVec, Mx{idxPulse});
                    case 'contour'
                        contour(posVec, csVec, Mx{idxPulse});
                end
                set(ax(counter), 'CLim', [-1 1]);
                title('M_x');
                response{idxPulse, idxPlotType} = Mx{idxPulse};
            case 'my'
                switch lower(plotType)
                    case '2d'
                        imagesc(posVec, csVec, My{idxPulse});
                    case 'mesh'
                        mesh(posVec, csVec, My{idxPulse});
                    case 'contour'
                        contour(posVec, csVec, My{idxPulse});
                end
                set(ax(counter), 'CLim', [-1 1]);
                title('M_y');
                response{idxPulse, idxPlotType} = My{idxPulse};
            case 'mxy'
                switch lower(plotType)
                    case '2d'
                        imagesc(posVec, csVec, abs(Mxy{idxPulse}));
                    case 'mesh'
                        mesh(posVec, csVec, abs(Mxy{idxPulse}));
                    case 'contour'
                        contour(posVec, csVec, abs(Mxy{idxPulse}));
                end
                set(ax(counter), 'CLim', [0 1]);
                title('|M_{xy}|');
                response{idxPulse, idxPlotType} = abs(Mxy{idxPulse});
            case 'atan'
                switch lower(plotType)
                    case '2d'
                        imagesc(posVec, csVec, atan2(My{idxPulse}, Mx{idxPulse}));
                    case 'mesh'
                        mesh(posVec, csVec, atan2(My{idxPulse}, Mx{idxPulse}));
                    case 'contour'
                        contour(posVec, csVec, atan2(My{idxPulse}, Mx{idxPulse}));
                end
                set(ax(counter), 'CLim', [-pi pi]);
                title('atan');
                response{idxPulse, idxPlotType} = atan2(My{idxPulse}, Mx{idxPulse});
            case 'atandeg'
                switch lower(plotType)
                    case '2d'
                        imagesc(posVec, csVec, atan2(My{idxPulse}, Mx{idxPulse})/pi*180);
                    case 'mesh'
                        mesh(posVec, csVec, atan2(My{idxPulse}, Mx{idxPulse})/pi*180);
                    case 'contour'
                        contour(posVec, csVec, atan2(My{idxPulse}, Mx{idxPulse})/pi*180);
                end
                set(ax(counter), 'CLim', [-180 180]);
                title('atan (deg.)');
                response{idxPulse, idxPlotType} = atan2(My{idxPulse}, Mx{idxPulse})/pi*180;
        end
        xlabel('mm');
        ylabel('kHz');
        colorbar
    end
end
switch lower(plotType)
    case '2d'
        colormap('gray');
    case 'mesh'
        % Natural color map
    case 'contour'
        % Natural color map
end
linkaxes(ax, 'x');

if numPulses==1
    set(gcf, 'Name', sprintf('Duration: %.2f (ms), Peak B1: %.2f (kHz), Num Steps: %d\n', pulses{1}.tp, max(pulses{1}.RFamp), numel(pulses{1}.RFamp)));
end

