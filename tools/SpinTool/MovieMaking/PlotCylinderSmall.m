function PlotCylinderSmall
% Plots just a white cylinder.

zi = -5;
zf = 5;

% Create and resize figure
figHnd = figure;
set(figHnd,'Position',[100 100 130 260]);
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
plot3(Xt,Yt,Zt,'k', 'LineWidth',1);
tAxis = linspace(-pi/2,pi/2,Nt);
Xt = cos(tAxis);
Yt = sin(tAxis);
Zt = ones(1,Nt)*zi;
plot3(Xt,Yt,Zt,'k', 'LineWidth',1);
tAxis = linspace(pi/2,3*pi/2,Nt);
Xt = cos(tAxis);
Yt = sin(tAxis);
Zt = ones(1,Nt)*zi;
plot3(Xt,Yt,Zt,'k--','LineWidth',1);
% Connect them with a line
tAxis = linspace(zi,zf,Nt);
Xt = zeros(1,Nt);
Yt = ones(1,Nt);
Zt = tAxis;
plot3(Xt,Yt,Zt,'k', 'LineWidth',1);
Xt = zeros(1,Nt);
Yt = -ones(1,Nt);
Zt = tAxis;
plot3(Xt,Yt,Zt,'k', 'LineWidth',1);



% Plot z-axis
axisHnd = plot3([0 0], [0 0], [zi-3 zf+3],'k-','LineWidth',1.5);

