function figHnd = PlotSpinsInCylinder(phaseVec, zAxis, vecColor, arrowThickness)
% SYNTAX: PlotSpinsInCylinder(phaseVec, zAxis, arrowThickness)
% Plots isochromats along a cylinder with given phases. For example, useful
% for plotting the spins wind and unwind during acquisition.
%   phaseVec - input vector of phases
%   zAxis - input vector of z-axis. Must have same number of elements as
%           phaseVec. 
%   vecColor - 3 element vector denoting the color of the arrows
%   arrowThickness - can be omitted. Thickness of arrows. Default is 0.015

if nargin<4
    arrowThickness = 0.015;
end;

% Number of slices. Should also equal the number of entries in phaseVec.
Nz = length(zAxis);

% First, create a 3xNz matrix, to hold the vectors
for k=1:Nz
    mag(k,:) = [cos(phaseVec(k)), sin(phaseVec(k)), 0];
end


% Cylinder length
zf = max(zAxis);
zi = min(zAxis);
L = abs(zf-zi);


% Create and resize figure
figHnd = figure;
set(figHnd,'Position',[100 100 200 400]);
figPos = get(figHnd,'Position');


% Plot Cylinder
ttt = 0:pi/10:2*pi;
[x,y,z] = cylinder(ones(1,50));
% cylHnd=surf(x,y,z*L+zi);
xlim=1.3;
r=1;
axis([-r*xlim r*xlim -r*xlim r*xlim zi zf]);
% set(cylHnd,'FaceAlpha',0);
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


% Adjust axis and figure properties to make background white
set(gcf,'color','white');
set(gca,'Visible','off');
campos([26.5730  0   13.9263]);

% Draw the outline of the cynlinder
% First, draw circles at top & bottom of cylinder
Nt = 100;
tAxis = linspace(0,2*pi,Nt);
Xt = cos(tAxis);
Yt = sin(tAxis);
Zt = ones(1,Nt)*zf;
plot3(Xt,Yt,Zt,'k', 'LineWidth',1.5);
tAxis = linspace(-pi/2,pi/2,Nt);
Xt = cos(tAxis);
Yt = sin(tAxis);
Zt = ones(1,Nt)*zi;
plot3(Xt,Yt,Zt,'k', 'LineWidth',1.5);
tAxis = linspace(pi/2,3*pi/2,Nt);
Xt = cos(tAxis);
Yt = sin(tAxis);
Zt = ones(1,Nt)*zi;
plot3(Xt,Yt,Zt,'k--','LineWidth',1.2);
% Connect them with a line
tAxis = linspace(zi,zf,Nt);
Xt = zeros(1,Nt);
Yt = ones(1,Nt);
Zt = tAxis;
plot3(Xt,Yt,Zt,'k', 'LineWidth',1.5);
Xt = zeros(1,Nt);
Yt = -ones(1,Nt);
Zt = tAxis;
plot3(Xt,Yt,Zt,'k', 'LineWidth',1.5);


% Compute positions of spins along the sample
Dz = L/Nz;
spinScale = 1;

% Plot z-axis
axisHnd = plot3([0 0], [0 0], [zi-2.5 zf+2.5],'k-','LineWidth',1.5);

% Plot spins ("arrows")
curVec = [0 0 0];
for p=1:Nz
    previousVec = curVec;
    curVec = [0 0 zAxis(p)] + [mag(p,1) mag(p,2) 0];
    hnd(p) = arrowPlot([0 0 zAxis(p)], curVec,vecColor, arrowThickness);
    if p~=1
        plot3([previousVec(1) curVec(1)],  [previousVec(2) curVec(2)], [previousVec(3) curVec(3)], 'color', vecColor);
    end
end;