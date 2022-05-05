function [movieSpins, spins, fid]=MovieBlochIsochromats(pulse,spins,numFrames,numPlottedSpins, plotRF, isCloseWindow, isInvertedColor, maxMovieDuration)
% SYNTAX: 
%
% [movieSpins, spins, fid] = MovieBlochIsochromats(pulse,
%                                                  spins,
%                                                  numFrames,
%                                                  numPlottedSpins, 
%                                                  plotRF)
%
% Creates a Matlab movie object with several isochromats, using the data
% in the input spin structure. Plots the Bloch sphere!
%
% Input Variable      Type     Range     Description
% pulse               pulse    N/A       RF pulse to be applied
% spins               spins    N/A       Spin structure to be visualized
% numFrames           double   1 - inf   Number of frames in movie. 
%                                        If > number of steps in pulse, 
%                                        the function will be aborted.
% numPlottedSpins     double   1 - inf   Number of spins to be plotted. 
%                                        If > number of spins, the function 
%                                        will be aborted.
% plotRF              boolean  0, 1      If set to true, the RF (red) will
%                                        be plotted.
% isCloseWindow       boolean  0, 1      If set to true, the movie window
%                                        will close automatically when
%                                        the movie ends
% isInvertedColor     boolean  0, 1      If set to true, the background
%                                        will be black
% maxMovieDuration    ms       0 to      The movie will be played back from
%                              pulse.tp  t=0 to t=maxMovieDuration.
%
% Note: you should keep the number of spins under 10, otherwise the movie
% will seem extremely cluttered.


% Vector containing the plotted spins' numberings. For example, if the same
% has 100 spins, and 10 are plotted, this vector will have the form:
%    1    12    23    34    45    56    67    78    89   100
% It may be non equi-spaced due to the constraint first and last spins 
% should be included.
numSpins = numel(spins);
if numPlottedSpins > numSpins
    disp('Too many plotted spins requested. Aborting ... ');
    disp(['(',num2str(numPlottedSpins),' requested, ', num2str(numSpins),' in sample.)']);
    return
end
plottedSpinsVec = floor(linspace(1,numSpins, numPlottedSpins));

% Find the maximal frame, as dictated by maxMovieDuration
maxFrameIdx = round(min(maxMovieDuration/pulse.tp, 1)*length(pulse.RFamp));

% Create vector 1 .. num of frames. It may be 
numTimeSteps = maxFrameIdx;
if numFrames > numTimeSteps 
    disp('Too many frames requested. Aborting ... ');
    disp(['(',num2str(numFrames),' requested, ', num2str(numTimeSteps),' steps in pulse.)']);
    return
end
frameVec = floor(linspace(1, numTimeSteps, numFrames));

spinColor = [0 0 1];   % Blue
RFColor   = [1 0 0];   % Red
% Assign different colors to the different plotted spins.
for curPlottedSpin=1:numPlottedSpins
    if isInvertedColor
        spinColorVec{curPlottedSpin} = spinColor +  (1-spinColor)*(numPlottedSpins - curPlottedSpin)/numPlottedSpins;
    else
        spinColorVec{curPlottedSpin} = spinColor - spinColor*(numPlottedSpins - curPlottedSpin)/numPlottedSpins;
    end
end

% Compute the evolution of the magnetization under the application of the
% given pulse. Extract relevant evolution.
for curPlottedSpin=1:numPlottedSpins
    idxCurrentSpin = plottedSpinsVec(curPlottedSpin);
    [Mx, My, Mz] = ApplyPulseRelaxDiagnostics(spins(idxCurrentSpin).cs, spins(idxCurrentSpin).r, spins(idxCurrentSpin).M, pulse, spins(idxCurrentSpin).T1, spins(idxCurrentSpin).T2, spins(idxCurrentSpin).M0);
    magTimeEvol{curPlottedSpin}.Mx = Mx(frameVec);
    magTimeEvol{curPlottedSpin}.My = My(frameVec);
    magTimeEvol{curPlottedSpin}.Mz = Mz(frameVec);
end;

% Compute B1(t), in kHz - just the RF (without the z component)
rf     = [pulse.RFamp.*cos(pulse.RFphase); pulse.RFamp.*sin(pulse.RFphase); 0.*zeros(1,length(pulse.RFamp))];
rfMax = max(max(abs(rf)));
rf    = rf./rfMax;  % Normalize to [0 1]

% Plot the 3D Bloch Sphere
hh=figure;
PlotBloch(isInvertedColor);
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

% Create movie
for curFrame = 1:numFrames
    % Plot spins
    for curPlottedSpin=1:numPlottedSpins
        M = [magTimeEvol{curPlottedSpin}.Mx(curFrame), ... 
             magTimeEvol{curPlottedSpin}.My(curFrame), ... 
             magTimeEvol{curPlottedSpin}.Mz(curFrame)];
        q1(curPlottedSpin)=arrowPlot([0 0 0],M,spinColorVec{curPlottedSpin},0.015);
    end;
    % Plot RF
    if plotRF == true
        q2=arrowPlot([0 0 0],rf(:,curFrame)',RFColor,0.015);
    end;
    movieSpins(curFrame) = getframe(hh,[1 1 QQ(3) QQ(4)]);
    % Delete plotted objects to free up memory
    if ((curFrame~=numFrames) || (isCloseWindow))
        for curPlottedSpin=1:numPlottedSpins
            delete(q1(curPlottedSpin));
            clear q1(curPlottedSpin)
        end;
        if plotRF == true
            delete(q2);
            clear q2
        end;
    end
end
if isCloseWindow, delete(hh); end

% Propagate spins and retrieve FID
[spins, fid] = AcquirePulseRelax(spins, pulse);
