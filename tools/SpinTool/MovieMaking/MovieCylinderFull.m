function [F,acqSignal, Mfinal]=MovieCylinderFull(pulse, NSpins, zi, zf, cs, Mi, NSteps, arrowThickness)
% Notes:
%   zf>zi

if nargin<8
    arrowThickness = 0.015;
end;


% If Mi is a single vector, "repmat" it so it becomes an 3xNSpins matrix
if size(Mi,2) == 1
    Mi = repmat(Mi,1,NSpins);
end;
% Compute maximal gradient - for colormap adjustments
GzMax = max(abs(pulse.Gz));
% Cylinder length
L = abs(zf-zi);
% Number of time steps in RF pulse
numTimeSteps = length(pulse.RFamp);


% Create and resize figure
figHnd = figure;
set(figHnd,'Position',[100 100 200 400]);
figPos = get(figHnd,'Position');


% Plot Cylinder
ttt = 0:pi/10:2*pi;
[x,y,z] = cylinder(ones(1,50));
cylHnd=surf(x,y,z*L+zi);
xlim=1.3;
r=1;
axis([-r*xlim r*xlim -r*xlim r*xlim zi zf]);
set(cylHnd,'FaceAlpha',0.4);
shading flat;
hold
grid off
load gradMap;

% Create Color maps
colorPos  = [1 1 0];   % Cyan
colorZero = [0 0.7 0.7];   % Yellow
colorNeg  = [1 0 0];   % Red
nColors = 32;
cylMapPos  = [linspace(colorNeg(1),colorZero(1),nColors), linspace(colorZero(1),colorPos(1),nColors);
              linspace(colorNeg(2),colorZero(2),nColors), linspace(colorZero(2),colorPos(2),nColors);
              linspace(colorNeg(3),colorZero(3),nColors), linspace(colorZero(3),colorPos(3),nColors)]';
cylMapNeg  = cylMapPos(end:-1:1,:);
cylMapZero = repmat(colorZero,nColors*2,1);
%set(gca,'CLim',[0.3 0.7]);


% Adjust axis and figure properties to make background white
set(gcf,'color','white');
set(gca,'Visible','off');
campos([26.5730  0   13.9263]);


% Compute positions of spins along the sample
Dz = L/(NSpins);
zz = linspace(zi, zf, NSpins);
%spinScale = Dz/1.5/1.1;
spinScale = 1;


% Compute evolutions of spins at corresponding positions
for p=1:NSpins
    [mag{p}.Mx, mag{p}.My, mag{p}.Mz] = ApplyPulseDiagnostics(cs, [0; 0; zz(p)], Mi(:,p), pulse);
    mag{p}.Mx = mag{p}.Mx.*spinScale;
    mag{p}.My = mag{p}.My.*spinScale;
    mag{p}.Mz = mag{p}.Mz.*spinScale;
end;

% Create movie
for k=1:NSteps:numTimeSteps
    % Compute acquired signal    
    for p=1:NSpins
        mxy(p) = mag{p}.Mx(k) + i*mag{p}.My(k);
    end;
    acqSignal(k) = sum(mxy);
    gradParam = abs(pulse.Gz(k))/GzMax;
    if pulse.Gz(k)>0
        set(gcf,'colormap',gradParam.*cylMapPos + (1-gradParam).*cylMapZero);
    else
        set(gcf,'colormap',gradParam.*cylMapNeg + (1-gradParam).*cylMapZero);
    end;
    % Plot z-axis
    axisHnd = plot3([0 0], [0 0], [zi-2.5 zf+2.5],'r-','LineWidth',2);
    % Plot spins ("arrows")
    for p=1:NSpins
        hnd(p) = arrowPlot([0 0 zz(p)], [0 0 zz(p)] + [mag{p}.Mx(k) mag{p}.My(k) mag{p}.Mz(k)],[0 0 0], arrowThickness);
    end;
    % Capture frame
    F(k) = getframe(figHnd,[1 1 figPos(3) figPos(4)]);
    % Remove spins
    for p=1:NSpins
        delete(hnd(p));
    end;
    delete(axisHnd);
end;


% Compute output magnetizations (Mfinal) after the pulse
for p=1:NSpins
    Mfinal(1,p) = mag{p}.Mx(end);
    Mfinal(2,p) = mag{p}.My(end);
    Mfinal(3,p) = mag{p}.Mz(end);
end;


