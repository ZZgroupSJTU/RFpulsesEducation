function outputMovie=MovieBlochB(B, pulseDuration, spinPos,initialMag,numFramesToSkip,addAnnotations, isNormGlobalB, isShowNormalPlane, isPlotSpinTrajectory, isInvertedColors, isPlotBlochSphere, isPlotProgressBar)
% Creates an animated movie of an RF pulse and a spin having a defined
% position and chemical shift.
%
% Inputs
%
% Variable Name     Units      Description
% spinPos           mm         Position of spin
% initialMag        -          Initial magnetization of spin. Can be 
%                              either normalized or unnormalized.
% numFramesToSkip   -          Number of pulse steps to skip between movie
%                              frames. Example: if a pulse has 150 steps,
%                              and numFramesToSkip = 1, the movie will have
%                              150 frames, which might be too long. 
%                              Selecting numFramesToSkip = 10 will only
%                              add to the movie steps 1, 11, 21, ..., 141]
%                              which might make more sense, visually and
%                              memory-wise.
% addAnnotation     -          True or false. If set to true, the position
%                              and chemical shift information for the 
%                              current spin will be added to each frame.
% isNormGlobalB     -          Since B(t) has units and M(t) doesn't,
%                              the movie will display a normalized B(t).
%                              If set to 0, B(t) will be normalized at each 
%                              time step, resulting in unity magnitude 
%                              throughout. If set to any other number,
%                              B(t) will be normalized to its global value
%                              times the number (e.g., set to 0.5 to 
%                              normalize B(t) to its global value, and then
%                              multiply each element by 0.5).
% isShowNormalPlane -          If set to 1, the plane normal to B(t) will
%                              be plotted as a function of time. If set to
%                              0, the plane will not be plotted.
% isPlotSpinTrajectory -       If set to 1, a "trajectory" (line) will be
%                              traced by the spin's end point will be 
%                              shown.
% isInvertedColors  -          If set to 1, the Bloch sphere will be bright
%                              against a black background.
% isPlotBlochSphere            Default: true. Set to false to eliminate
%                              Bloch sphere.
% isPlotProgressBar            Default: true.

if nargin<11, isPlotBlochSphere = true; end
if nargin<12, isPlotProgressBar = true; end

% ------------------------------------------------------------------------
% Verify inputs 
% ------------------------------------------------------------------------

% Normalize input magnetization

numPulseSteps = size(B,2);
if (isInvertedColors)
    spinColor = [0.4 0.4 1];   % Blue
    RFColor   = [1 0.2 0.2];   % Red
else
    spinColor = [0 0 1];   % Blue
    RFColor   = [1 0 0];   % Red
end


% Vector specifying which time steps will make it into the final movie
frameVector = [1:numFramesToSkip:numPulseSteps];
if (frameVector(end)~=numPulseSteps)
    frameVector(end+1) = numPulseSteps;
end
numFrames = length(frameVector);


% ------------------------------------------------------------------------
% Simulate spin's evolution
% ------------------------------------------------------------------------

M = initialMag;
dwellTime = pulseDuration/numPulseSteps;
for idx=1:numPulseSteps
    Mx(idx) = M(1);
    My(idx) = M(2);
    Mz(idx) = M(3);
    Bt = B(:,idx);
    rotAngle = 2*pi*norm(Bt)*dwellTime;
    M = RotMat(Bt, rotAngle)*M;
end

Mt=[Mx(frameVector); My(frameVector); Mz(frameVector)];


% ------------------------------------------------------------------------
% Normalize RF vector to a unit (direction) vector
% ------------------------------------------------------------------------

if isNormGlobalB>0
    B = B./max(sqrt(sum(B.*B,1))).*isNormGlobalB;  
else
    B = B./repmat(sqrt(sum(B.*B,1)), 3,1);  
end

% Takes only those values of B(t) that will be featured in the movie
Bft = [B(1, frameVector); B(2, frameVector); B(3,frameVector)];

% ------------------------------------------------------------------------
% Set the stage: draw the 3D Bloch Sphere and progress bars
% ------------------------------------------------------------------------

% Initialize main movie figure, and plot the Bloch sphere
mainFigHandle=figure;
if isPlotBlochSphere
    PlotBloch(isInvertedColors);
else
    PlotBloch(isInvertedColors, 0, [], 0);
end
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
if isInvertedColors
    set(gcf,'color','black');
    set(gca,'color','black');
else
    set(gcf,'color','white');
    set(gca,'color','white');
end
set(gca,'LineStyleOrder','--');
set(gca,'Visible','off');
figPos = get(mainFigHandle,'Position');

% Draw progress bar
if isPlotProgressBar
    progressBarWidth = 0.6;
    progressBarPosX=0.22;
    progressBarPosY=0.12;
    annotation('rectangle',[progressBarPosX progressBarPosY-0.01 progressBarWidth 0.02],'EdgeColor','blue','LineWidth',1);
    str2(1) = {'Progress Bar'};
    str2(2) = {'t=0'};
    str2(3) = {['t=',num2str(pulseDuration),' ms']};
end
    
% Add annotations: write auxilary info (chemShift, progress, position and so forth)
if addAnnotations == true
    text(5.75,0.8,str2(1),'HorizontalAlignment','right')
    text(7.1,-0.15,str2(2),'HorizontalAlignment','right')
    text(6.9,1.8,str2(3),'HorizontalAlignment','right')
    text(-6,-1.81,'Position:','HorizontalAlignment','right')
    text(-5,-1.6,['x: ',num2str(spinPos(1)),' mm'],'HorizontalAlignment','right')
    text(-4,-1.5,['y: ',num2str(spinPos(2)),' mm'],'HorizontalAlignment','right')
    text(-3,-1.4,['z: ',num2str(spinPos(3)),' mm'],'HorizontalAlignment','right')
end;


% ------------------------------------------------------------------------
% Play out movie and capture the screen's output
% ------------------------------------------------------------------------

% Initialize output for speed
outputMovie(numFrames) = getframe(mainFigHandle,[1 1 figPos(3) figPos(4)]);

% Create movie
for curFrame = 1:numFrames
    if isShowNormalPlane
        q4=DrawCircle3D(gca, 0.5, B(:,curFrame));
    end
    q1=arrowPlot([0 0 0],Mt(:,curFrame)',spinColor,0.015);
    q2=arrowPlot([0 0 0],Bft(:,curFrame)',RFColor,0.015);
    if isPlotProgressBar
        q3=annotation('line', [progressBarPosX progressBarPosX+progressBarWidth*(curFrame-1)/(numFrames-1)], [progressBarPosY progressBarPosY],'LineWidth',5,'Color','red');
    end
    if isPlotSpinTrajectory
        q5 = line(Mt(1,1:curFrame), Mt(2,1:curFrame), Mt(3,1:curFrame), 'Color', spinColor);
    end
    
    
    pause(0.01);
    outputMovie(curFrame) = getframe(mainFigHandle,[1 1 figPos(3) figPos(4)]);

    if (curFrame~=numFrames)
        delete(q1);
        delete(q2);
        if isPlotProgressBar
            delete(q3);
        end
        if isShowNormalPlane
            if (q4~=0)
                delete(q4);
            end
        end
        if (isPlotSpinTrajectory)
            if q5~=0
                delete(q5)
            end
        end
    end
end

