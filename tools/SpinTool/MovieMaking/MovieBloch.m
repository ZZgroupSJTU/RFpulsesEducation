function outputMovie=MovieBloch(pulse,pulseParameter,spinPos,initialMag,numFramesToSkip,addAnnotations, plotJustRF, isPlotPlaneRF, cameraPosition)
% SYNTAX:
%
% outputMovie=MovieBloch(pulse, pulseParameter, spinPos, initialMag,
%                        numFramesToSkip, addAnnotations, plotJustRF, 
%                        isPlotPlaneRF, cameraPosition)
%
% Creates an animated movie of an RF pulse and a spin having a defined
% position and chemical shift.
% HINT: 
%
% Inputs
%
% Variable Name     Units      Description
% pulse             -          Input pulse structure 
%                              OR
%                              a 3xN matrix containing B(t)
% pulseParameter    kHz        Chemical shift of spin (if pulse is a struct)
%                              OR
%                   ms         Duration of RF field (if pulse is a matrix)
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
% plotJustRF        -          If set to true, just the RF vector will
%                              be plotted, and the offset omitted. Setting
%                              this to false plots the entire B field in
%                              every frame.
% isPlotPlaneRF     -          True/false. If true, plots a plane (circle)
%                              perpendicular to the instantaneous field.
% cameraPosition    -          1x3 vector denoting camera view


% ------------------------------------------------------------------------
% Verify inputs 
% ------------------------------------------------------------------------

if (isstruct(pulse))
    numPulseSteps = length(pulse.RFamp);
    pulseDuration = pulse.tp;
    chemShift = pulseParameter;
else
    pulseDuration = pulseParameter;
    numPulseSteps = size(pulse, 2);
    chemShift = 0;
end

% Set colors of B(t) and M(t) vectors
spinColor = [0 0 1];   % Blue
RFColor   = [1 0 0];   % Red


% Vector specifying which time steps will make it into the final movie
frameVector = [1:numFramesToSkip:numPulseSteps];
if (frameVector(end)~=numPulseSteps)
    frameVector(end+1) = numPulseSteps;
end
numFrames = length(frameVector);



% ------------------------------------------------------------------------
% Simulate spin's evolution
% ------------------------------------------------------------------------

Mx = zeros(1, numPulseSteps);
My = zeros(1, numPulseSteps);
Mz = zeros(1, numPulseSteps);

if (isstruct(pulse))
    [Mx, My, Mz] = ApplyPulseDiagnostics(chemShift, spinPos, initialMag, pulse);
else
    M = initialMag;
    dwellTime = pulseDuration/numPulseSteps;
    for idx=1:numPulseSteps
        Mx(idx) = M(1);
        My(idx) = M(2);
        Mz(idx) = M(3);
        Bt = pulse(:,idx);
        rotAngle = 2*pi*norm(Bt)*dwellTime;
        M = RotMat(Bt, rotAngle)*M;
    end
end

Mt=[Mx(frameVector); My(frameVector); Mz(frameVector)];


% ------------------------------------------------------------------------
% Compute RF vector
% ------------------------------------------------------------------------

if (isstruct(pulse))
    Bx = pulse.RFamp(frameVector).*cos(pulse.RFphase(frameVector));
    By = pulse.RFamp(frameVector).*sin(pulse.RFphase(frameVector));
    Bz = pulseParameter.*ones(1,numFrames);
else
    Bx = pulse(1,frameVector);
    By = pulse(2,frameVector);
    Bz = pulse(3,frameVector);
end



if (plotJustRF)
    B = [Bx; By; 0.*Bz];
else
    B = [Bx; By; Bz];
end

B = B./max(sqrt(sum(B.*B,1)));



% ------------------------------------------------------------------------
% Set the stage: draw the 3D Bloch Sphere and progress bars
% ------------------------------------------------------------------------

% Initialize main movie figure, and plot the Bloch sphere
mainFigHandle=figure;
PlotBloch;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
set(gcf,'color','white');
set(gca,'color','white');
set(gca,'LineStyleOrder','--');
set(gca,'Visible','off');
if (nargin>=9)
    set(gca,'CameraPosition', cameraPosition);
end
figPos = get(mainFigHandle,'Position');

% Draw progress bar
progressBarWidth = 0.6;
progressBarPosX=0.22;
progressBarPosY=0.12;
annotation('rectangle',[progressBarPosX progressBarPosY-0.01 progressBarWidth 0.02],'EdgeColor','blue','LineWidth',1);
str2(1) = {'Progress Bar'};
str2(2) = {'t=0'};
str2(3) = {['t=',num2str(pulseDuration),' ms']};

% Add annotations: write auxilary info (chemShift, progress, position and so forth)
if addAnnotations == true
    text(5.75,0.8,str2(1),'HorizontalAlignment','right')
    text(7.1,-0.15,str2(2),'HorizontalAlignment','right')
    text(6.9,1.8,str2(3),'HorizontalAlignment','right')
    if (isstruct(pulse))
        text(-7,-1.78,['Offset: ',num2str(pulseParameter),' kHz'],'HorizontalAlignment','right')
    end
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
    if (isPlotPlaneRF)
        q4=DrawCircle3D(gca, 0.5, B(:,curFrame));
    end
    q1=arrowPlot([0 0 0],Mt(:,curFrame)',spinColor,0.015);
    q2=arrowPlot([0 0 0],B(:,curFrame)',RFColor,0.015);
    q3=annotation('line', [progressBarPosX progressBarPosX+progressBarWidth*(curFrame-1)/(numFrames-1)], [progressBarPosY progressBarPosY],'LineWidth',5,'Color','red');
    
    pause(0.01);
    outputMovie(curFrame) = getframe(mainFigHandle,[1 1 figPos(3) figPos(4)]);

    if (curFrame~=numFrames)
        delete(q1);
        delete(q2);
        delete(q3);
        if (isPlotPlaneRF)
            if (q4~=0)
                delete(q4);
            end
        end
    end
end

