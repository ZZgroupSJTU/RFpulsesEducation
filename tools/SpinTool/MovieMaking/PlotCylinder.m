function figHandle=PlotCylinder(height,radius,facth,factr);
% Plots a 3D Cylinder, and returns the handle to the figure.

figHandle = figure;

% Create axes
hold on

% set axes
axis([-radius*factr radius*factr -radius*factr radius*factr -height/2*facth height/2*facth]);

% Create "sphere"
N=50;
x=linspace(-1,1,N);
y1 = sqrt(1-x.^2);
y2 = -y1;

plot3(x,y1,height/2.*ones(1,N),'k','LineWidth',1);
plot3(x,y2,height/2.*ones(1,N),'k','LineWidth',1);

plot3(x,y1,-height/2.*ones(1,N),'k','LineWidth',1);
plot3(x,y2,-height/2.*ones(1,N),'k','LineWidth',1);

plot3([radius radius], [0 0], [-height/2, height/2], 'k', 'LineWidth',1);
plot3([-radius -radius], [0 0], [-height/2, height/2], 'k', 'LineWidth',1);


campos([0.99 0.1 0.2])