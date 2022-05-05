function [movieSpins, movieTimeAxis, fid]=MovieBlochIsochromatsSeq(seq, varargin)
% function [movieSpins, movieTimeAxis, fid]=MovieBlochIsochromatsSeq(seq,spins,numFrames,numPlottedSpins, isPlotRF, isCloseWindow, isInvertedColor, isPlotFID, isPlotSpinTrajectories, maxMovieDuration, spinColorVec)
% 
% [movieSpins, movieTimeAxis, fid]=MovieBlochIsochromatsSeq(seq, ...)
% 
% Creates a Matlab movie object with several isochromats, using the data
% in the input spin structure. Plots the Bloch sphere!
%
% Optional argument-value pairs:
%
%   spins                    Spin structure to be visualized
%   numFrames                Number of frames in movie. If > number of 
%                            steps in pulse, the function will be aborted.
%   numPlottedSpins          Number of spins to be plotted. If > number of 
%                            spins, the function will be aborted.
%   isPlotRF                 If set to true, the RF (red) will be plotted.
%   isPlotSpinsXY            Adds a plot of the spins in the xy-plane
%                            (a "bird's eye view" - easier to visualize phases)
%   isCloseWindow            If set to true, the movie window will close 
%                            automatically when the movie ends.
%   isInvertedColor          If set to true, the background will be black
%   isPlotFID                If set to 1, ALL of the spins will be 
%                            simulated and the FID will be calculated 
%                            (Mx+1i*My) at each time point throughout the 
%                            evolution & plotted on the screen.
%   isPlotMz                 
%   title                    An optional string title
%   isPlotSpinTrajectories   If set to 1, trajectories of each isochromat 
%                            will be drawn on the surface of the Bloch
%                            sphere.
%   maxMovieDuration         The movie will be played back from t=0 to 
%                            t=maxMovieDuration (in ms).
%   spinColorVec             Cell array of 1x3 vectors of RGB color values 
%                            for each spin in the isochromat. If omitted, 
%                            an automatic matrix will be generated.
%
% Note: you should keep the number of spins under ~10, otherwise the movie
% will seem extremely cluttered.

% ========================================================================
% Parse Input
% ========================================================================

isBoolean = @(x)islogical(x) || isnumeric(x) && all(x(:)==0 | x(:)==1);

p = inputParser;
p.addRequired('seq', @(x) IsSequence(x));
p.addOptional('spins', InitSpinsRelax(0, 1, 1, [0;0;1], 1e6, 1e6, 1));
p.addOptional('numFrames', 0);
p.addOptional('numPlottedSpins', 1);
p.addOptional('isPlotRF', false);
p.addOptional('isPlotSpinsXY', false);
p.addOptional('isCloseWindow', false);
p.addOptional('isInvertedColor', false, @(x) isBoolean(x));
p.addOptional('isPlotFID', false, @(x) isBoolean(x));
p.addOptional('isPlotSpinTrajectories', false, @(x) isBoolean(x));
p.addOptional('maxMovieDuration', CalcSeqTotalTime(seq), @(x) isnumeric(x));
p.addOptional('isPlotMz', true, @(x) isBoolean(x));
p.addOptional('title', '', @(x) ischar(x));
p.addOptional('spinColorVec', []);

p.parse(seq, varargin{:});

spins = p.Results.spins;
numFrames = p.Results.numFrames;
numPlottedSpins = p.Results.numPlottedSpins;
isPlotRF = p.Results.isPlotRF;
isCloseWindow = p.Results.isCloseWindow;
isInvertedColor = p.Results.isInvertedColor;
isPlotFID = p.Results.isPlotFID;
isPlotMz = p.Results.isPlotMz;
isPlotSpinTrajectories = p.Results.isPlotSpinTrajectories;
isPlotSpinsXY = p.Results.isPlotSpinsXY;
maxMovieDuration = p.Results.maxMovieDuration;
spinColorVec = p.Results.spinColorVec; 
movieTitle = p.Results.title;

if numFrames==0, numFrames = 1; end

% ========================================================================
% Simulate time evolution under supplied sequence
% ========================================================================


% Vector containing the plotted spins' numberings. For example, if the same
% has 100 spins, and 10 are plotted, this vector will have the form:
%    1    12    23    34    45    56    67    78    89   100
% It may be non equi-spaced due to the constraint first and last spins 
% should be included.
numSpins = numel(spins);
if numPlottedSpins > numSpins
    disp('Too many plotted spins requested.');
    disp(['(',num2str(numPlottedSpins),' requested, ', num2str(numSpins),' in sample.)']);
    fprintf('Setting number of plotted spins to %d\n', numSpins);
    numPlottedSpins = numSpins;
end
plottedSpinsVec = floor(linspace(1,numSpins, numPlottedSpins));

% Compute the evolution of the magnetization under the application of the
% given sequence. Extract relevant evolution.
numReps = 1; % The system will not be driven into equilibrium by repeating the sequence multiple times
if isPlotFID || isPlotMz % Simulate ALL spins
    for idxCurrentSpin=1:numSpins
        [Mx{idxCurrentSpin}, My{idxCurrentSpin}, Mz{idxCurrentSpin}, timeAxis, B1Amp, B1Phase] = ApplySequenceDiagnostics(seq, numReps, spins(idxCurrentSpin).cs, spins(idxCurrentSpin).r, spins(idxCurrentSpin).M, spins(idxCurrentSpin).T1, spins(idxCurrentSpin).T2, spins(idxCurrentSpin).M0);
    end;
else % Simulate only those spins that are going to be plotted
    for curPlottedSpin=1:numPlottedSpins
        idxCurrentSpin = plottedSpinsVec(curPlottedSpin);
        [Mx{curPlottedSpin}, My{curPlottedSpin}, Mz{curPlottedSpin}, timeAxis, B1Amp, B1Phase] = ApplySequenceDiagnostics(seq, numReps, spins(idxCurrentSpin).cs, spins(idxCurrentSpin).r, spins(idxCurrentSpin).M, spins(idxCurrentSpin).T1, spins(idxCurrentSpin).T2, spins(idxCurrentSpin).M0);
    end;
end

% Find the maximal frame, as dictated by maxMovieDuration
maxTimeIdx = find(timeAxis>=maxMovieDuration, 1, 'first');
if isempty(maxTimeIdx)
    maxTimeIdx = numel(timeAxis);
end
if numFrames > maxTimeIdx 
    disp('Too many frames requested. ');
    disp(['(',num2str(numFrames),' requested. There are only ', num2str(maxTimeIdx),' steps in requested time window.)']);
    fprintf('Setting number of frames to max.\n');
    numFrames = maxTimeIdx;
end

% Create vector 1 .. num of frames. 
frameVec = floor(linspace(1, maxTimeIdx, numFrames));


movieTimeAxis = timeAxis(frameVec);
fid = zeros(1, numFrames);
mz = zeros(1, numFrames);
if isPlotFID || isPlotMz
    for idxSpin=1:numSpins
        magTimeEvol{idxSpin}.Mx = Mx{idxSpin}(frameVec);
        magTimeEvol{idxSpin}.My = My{idxSpin}(frameVec);
        magTimeEvol{idxSpin}.Mz = Mz{idxSpin}(frameVec);
        fid = fid + magTimeEvol{idxSpin}.Mx + 1i*magTimeEvol{idxSpin}.My;
        mz = mz + magTimeEvol{idxSpin}.Mz;
    end
else
    for curPlottedSpin=1:numPlottedSpins
        magTimeEvol{curPlottedSpin}.Mx = Mx{curPlottedSpin}(frameVec);
        magTimeEvol{curPlottedSpin}.My = My{curPlottedSpin}(frameVec);
        magTimeEvol{curPlottedSpin}.Mz = Mz{curPlottedSpin}(frameVec);
        fid = fid + magTimeEvol{curPlottedSpin}.Mx + 1i*magTimeEvol{curPlottedSpin}.My;
        mz = mz + magTimeEvol{curPlottedSpin}.Mz;
    end;
end

spinColor = [0 0 1];   % Blue
RFColor   = [1 0 0];   % Red
% Assign different colors to the different plotted spins.
if isempty(spinColorVec)
    for curPlottedSpin=1:numPlottedSpins
        spinColorVec{curPlottedSpin} = spinColor - spinColor*(numPlottedSpins - curPlottedSpin)/numPlottedSpins;
    end
end

if isInvertedColor
    for curPlottedSpin=1:numPlottedSpins
        spinColorVec{curPlottedSpin} = [1 1 1] - spinColorVec{curPlottedSpin};
    end
end



% Plot the 3D Bloch Sphere
hFig=figure('Units', 'Pixels', 'Position', [300 300 800 600]);
axesBloch = axes('Parent', hFig);
if isPlotFID || isPlotMz
    axesFID = axes('Parent', hFig);
    axesFID.Units = 'Pixels';
    axesFID.XTick = [];
    axesFID.YTick = [];
    axesFID.ZTick = [];
    axesFID.LineStyleOrder = '--';
    axesFID.Visible = 'off';
    % set axes
    L = 1.1;
    axesFID.XLim = [-L L];
    axesFID.YLim = [-L L];
    axesFID.ZLim = [-L L];
    text(-1.05, 0.0, 0, 'RF', 'Parent', axesFID, 'FontSize', 14, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom');
    if isPlotFID && isPlotMz
        % RF + FID + Mz
    elseif isPlotFID
        % Just RF + FID
    else
        % Just RF + Mz
        text(-1.05, -1, 0, 'M_z', 'Parent', axesFID, 'FontSize', 14, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom');
        line([-1 1], [0 0], [0 0], 'LineStyle', '--', 'Color', 'black', 'Parent', axesFID);
    end
end
if isPlotSpinsXY
    axesPlane = axes('Parent', hFig);
    axesPlane.Units = 'Pixels';
    axesPlane.XTick = [];
    axesPlane.YTick = [];
    axesPlane.ZTick = [];
    axesPlane.LineStyleOrder = '--';
    axesPlane.Visible = 'off';
    % set axes
    L = 1.1;
    axesPlane.XLim = [-L L];
    axesPlane.YLim = [-L L];
    axesPlane.ZLim = [-L L];
end

[sphereCenter, scaleFactor] = PlotBloch(isInvertedColor, false, [], true, axesBloch);

set(axesBloch,'XTick',[]);
set(axesBloch,'YTick',[]);
set(axesBloch,'ZTick',[]);
if isInvertedColor
    set(hFig,'color','black');
    set(axesBloch,'color','black');
    if isPlotFID, axesFID.Color = 'black'; end 
    if isPlotSpinsXY, axesPlane.Color = 'black'; end
else
    set(hFig,'color','white');
    set(axesBloch,'color','white');
    if isPlotFID, axesFID.Color = 'white'; end 
    if isPlotSpinsXY, axesPlane.Color = 'white'; end
end
set(axesBloch,'LineStyleOrder','--');
set(axesBloch,'Visible','off');
axesBloch.Units = 'Pixels';

if ~isPlotSpinsXY && ~(isPlotFID || isPlotMz)
    % Just the Bloch sphere: Make the axes occupy the whole figure
    axesBloch.Position = [0 0 800 600];
    hText = text(0, 0, 1.3, movieTitle, 'HorizontalAlignment','center');    
elseif (isPlotSpinsXY && ~(isPlotFID || isPlotMz))
    % Bloch sphere + projection of spins onto XY plane
    axesBloch.Position = [0   150 400 300];
    axesPlane.Position = [400 150 400 300];
elseif (~isPlotSpinsXY && (isPlotFID || isPlotMz))
    % Bloch sphere + signal
    axesBloch.Position = [125 200 550 400];
    axesFID.Position = [50 50 700 150];
    hText = text(0, 0, 1.2, movieTitle, 'HorizontalAlignment','center', 'FontWeight', 'bold', 'Parent', axesBloch);    
else 
    % Plot both FID & projection of spins in the XY plane in addition to Bloch sphere
    axesBloch.Position = [0   250 400 300];
    axesPlane.Position = [400 250 400 300];
    axesFID.Position = [50 50 700 500];
    
end

if isPlotSpinsXY
    % Create "unit circle"
    axes(axesPlane);
    N=50;
    x=linspace(0,2*pi,N);
    line(cos(x), sin(x), 'Color',[1 1 1]*0.5,'LineWidth',1, 'LineStyle', ':');
end

QQ = get(hFig,'Position');



% Compute B1(t), in kHz - just the RF (without the z component)
rf    = [B1Amp.*cos(B1Phase); B1Amp.*sin(B1Phase); 0.*zeros(1,length(B1Phase))];
rfMax = max(max(abs(rf)));
rf    = rf./rfMax*scaleFactor; % Normalize to [0 scaleFactor]

scaledTimeAxis = (movieTimeAxis./max(movieTimeAxis) - 0.5)*2;
scaledFID = fid./max(abs(fid));
scaledRF = B1Amp(frameVec)./max(abs(B1Amp(frameVec)));
scaledMz = mz./max(abs(mz));

drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0 )  ;

% Create movie
for curFrame = 1:numFrames
    % Plot spins
    axes(axesBloch);
    for curPlottedSpin=1:numPlottedSpins
        if isPlotFID || isPlotMz
            idxCurSpin = plottedSpinsVec(curPlottedSpin);
        else
            idxCurSpin = curPlottedSpin; 
        end
        M = [magTimeEvol{idxCurSpin}.Mx(curFrame), ... 
             magTimeEvol{idxCurSpin}.My(curFrame), ... 
             magTimeEvol{idxCurSpin}.Mz(curFrame)];
        magM = norm(M); 
        M = M.*scaleFactor; 
        spinArrowSize = 0.015/magM;
        q1(curPlottedSpin)=arrowPlot(sphereCenter+[0 0 0],sphereCenter+M,spinColorVec{curPlottedSpin},spinArrowSize);
        if isPlotSpinTrajectories
            q5(curPlottedSpin) = line(sphereCenter(1)+magTimeEvol{curPlottedSpin}.Mx(1:curFrame)*scaleFactor, ...
                                      sphereCenter(2)+magTimeEvol{curPlottedSpin}.My(1:curFrame)*scaleFactor, ...
                                      sphereCenter(3)+magTimeEvol{curPlottedSpin}.Mz(1:curFrame)*scaleFactor, ...
                                      'Color', spinColorVec{curPlottedSpin}, 'LineStyle', '-');
        end
    end;
    % Plot RF
    if isPlotRF == true
        q2=arrowPlot(sphereCenter+[0 0 0],sphereCenter+rf(:,frameVec(curFrame))',RFColor,0.015);
    end;
    % Plot the FID
    if isPlotFID || isPlotMz
        axes(axesFID);
        q3RF = line(scaledTimeAxis(1:curFrame), real(scaledRF(1:curFrame)), 'LineStyle', '-', 'LineWidth', 2, 'Color', 'red'); 
        if isPlotFID && isPlotMz
            % Plot RF + FID + Mz
            q3Re = line(scaledTimeAxis(1:curFrame), real(scaledFID(1:curFrame)).*0.09 - 0.77, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'blue'); 
            q3Im = line(scaledTimeAxis(1:curFrame), imag(scaledFID(1:curFrame)).*0.09 - 1.0, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'magenta'); 
            q3Mz = line(scaledTimeAxis(1:curFrame), scaledMz(1:curFrame)*0.09 - 1.0, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black'); 
        elseif isPlotFID
            % Plot RF + FID (Real, Imag)
            q3Re = line(scaledTimeAxis(1:curFrame), real(scaledFID(1:curFrame)).*0.09 - 0.77, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'blue'); 
            q3Im = line(scaledTimeAxis(1:curFrame), imag(scaledFID(1:curFrame)).*0.09 - 1.0, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'magenta'); 
        elseif isPlotMz
            % Plot RF + Mz
            q3Mz = line(scaledTimeAxis(1:curFrame), scaledMz(1:curFrame)*0.2 - 0.85, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black'); 
        end
    end
    % Plot projections of spins in the xy-plane
    if isPlotSpinsXY
        axes(axesPlane)
        for curPlottedSpin=1:numPlottedSpins
            if isPlotFID
                idxCurSpin = plottedSpinsVec(curPlottedSpin);
            else
                idxCurSpin = curPlottedSpin; 
            end
            M = [magTimeEvol{idxCurSpin}.Mx(curFrame), ... 
                 magTimeEvol{idxCurSpin}.My(curFrame), ... 
                 0];
            magM = norm(M); 
            M = M.*scaleFactor; 
            spinArrowSize = 0.015/magM;
            q6(curPlottedSpin)=line([0 M(1)], [0 M(2)], 'LineStyle', '-', 'LineWidth', 2, 'Color', spinColorVec{curPlottedSpin}); 
        end;
        
    end
    
    % Plot isochromat trajectories
    movieSpins(curFrame) = getframe(hFig,[1 1 QQ(3) QQ(4)]);
    % Delete plotted objects to free up memory
    if ((curFrame~=numFrames) || (isCloseWindow))
        for curPlottedSpin=1:numPlottedSpins
            delete(q1(curPlottedSpin));
            clear q1(curPlottedSpin)
        end;
        if isPlotSpinTrajectories
            for curPlottedSpin=1:numPlottedSpins
                delete(q5(curPlottedSpin));
                clear q5(curPlottedSpin)
            end
        end
        if isPlotRF == true
            delete(q2);
            clear q2
        end;
        if isPlotSpinsXY
            for curPlottedSpin=1:numPlottedSpins
                delete(q6(curPlottedSpin));
                clear q6(curPlottedSpin)
            end;
        end
        if isPlotFID || isPlotMz
            delete(q3RF);
            clear q3RF;
            if isPlotFID
                delete(q3Re);
                delete(q3Im);
                clear q3Re;
                clear q3Im;
            end
            if isPlotMz
                delete(q3Mz);
                clear q3Mz;
            end
        end
    end
end
if isCloseWindow, delete(hFig); end

