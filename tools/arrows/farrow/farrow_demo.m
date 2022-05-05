% This script runs 5 'farrow' examples, in order

% Make sure to have 'farrow' function file in same folder

% Script:
% - Setup: Set figure and axes properties.
% - (1) Rainbow-Spiral-Arrows
% - (2) Arrows Along Circle
% - (3) Drawing Arrow
% - (4) Arrow with Cone-Patch
% - (5) Flower Arrow

% Setup

clear; clc; close all;

set(gcf, 'color', [.7 .7 .8], ...                % Figure, color
       'Menubar', 'none', ...
   'numbertitle', 'off');

set(gca,         'Units',  'Normalized', ...     % Axes  
         'OuterPosition', [0, 0, 1, .9], ...     % - position
  'TickLabelinterpreter',       'Latex', ...     % - font
  'XLim', [0 12], 'YLim',[-1 1], 'ZLim',[-1 1]); % - limits

info = annotation(  'TextBox',              ...  % Info Box
                   'Position', [0,.9,1,.1], ...  % - position
            'BackgroundColor',  [.2 .2 .2], ...  % - color
'Fontsize', 18, 'Interpreter',     'Latex');     % - font
                
view(120, 23);  grid on;                         % View Angle, 

% (1) Rainbow-Spiral-Arrows

set(info, 'string', 'Rainbow-Spiral Arrows',...
           'color', [.5 .5 .5] + rand(1, 3)/2);
       
for i = 0:.2:12                                 % Loop for each arrow
    r     = .5 + cos(i)*cos(i)/2;               % - 'color sphere'
    g     = .5 + cos(i)*sin(i)/2; 
    b     = .5 + sin(i)/2;
    color = [r, g, b];
    
    farrow(i, 0, 0, i, cos(i), sin(i), color);   % Arrow Function
    pause(.01);
end
pause(1);

% (2) Arrows along Circle

cla;

set(info, 'string', 'Arrows along Circle',...   % Update Info Text
     'color', [.5 .5 .5] + rand(1, 3)/2);       % New text color

set(gca, 'XLim', [-1 1], ...                    % New limits
         'YLim', [-1 1], ...
         'ZLim', [-1 1]);
     
line(cosd(0:360), sind(0:360), zeros(1, 361));  % Circle

for i = 0:5:360                                 % Loop around Circle:
    x = cosd(i);                                % - (x, y) follow circle
    y = sind(i);                              
    z = cosd(2*i);                              % - height up/down
    
    color = [0 0 1] + [.9 .9 0]*abs(z);         % - arrow color
    width = sqrt(10*abs(z));                    % - arrow width
    farrow(x, y, 0, x, y, z, color, width);     
    pause(.01);
end

pause(1);

% (3) Drawing Arrow

cla;                                            % Clear Axes
set(info, 'string','Drawing Arrow',...          % Update Info Text
    'color',   [.5 .5 .5] + rand(1, 3)/2);      % -- new text color
       
set(gca, 'XLim', [0 2], ...                     % New Axes Limits
         'YLim', [0 2], ...
         'ZLim', [0 2]);

line0 = animatedline('linewidth', 2);           % Init Animated Line

for i = 0:3:360                                 % Loop (upper curve)
    if i ~= 0, delete(a); end                   % - delete old arrow
    x = cosd(2*i + 90) + 1;                     % - define x,y,z
    y = i/180;
    z = 1.1;
    addpoints(line0, x, y, z);                  % - add to animated line
    a = farrow(.5, .5, .1, x, y, z);            % - new arrow
    pause(.01);
end

while z > 0                                     % Loop (straight down)
    z = z - .02;                                % - lower z
    if i ~= 0, delete(a); end                   % - delete old arrow
    
    addpoints(line0, x, y, z);                  % - add to animate line
    a = farrow(.5, .5, .1, x, y, z);            % - new arrow
    pause(.01);
end

for i = 0:2:360                                 % Loop (lower curve)
    delete(a); j = i + 180;                     % - delete old arrow
    
    x = cosd(j)*cosd(5*j + 90) + 1;             % - define (x, y)
    y = 2 - i/180;
    
    addpoints(line0, x, y, z);                  % - add to animated line
    a = farrow(.5, .5, .1, x, y, z);            % - new arrow
    pause(.01);
end

while z < 1.1                                   % Loop (down):
    delete(a);                                  % - delete arrow
    z = z + .03;                                % - move z down
    addpoints(line0, x, y, z);                  % - add to animated line
    a = farrow(.5, .5, .1, x, y, z);            % - new arrow 
    
    pause(.01);
end

pause(1); 

% (4) Arrow with Cone-Patch

cla;                                            % Clear Axes

set(info,'String', 'Arrow with Cone-Patch',...  % Update Text
         'Color', [.5 .5 .5] + rand(1, 3)/2);   % Update Text Color

set(gca, 'XLim', [-1 1], ...                    % New Axes Limits
         'YLim', [-1 1], ...
         'ZLim', [-1 1]);

light( 'Position',[1, 1, 0],'Style','local');   % Add Light   
     
faces    = [1:3:3000; 2:3:3000; 3:3:3000]';     % Cone Patch Face Data
vertices = ones(9000, 3);                       % Cone Patch Vertex Data

spiral = patch( ...                             % Cone Patch Object                 
        'faces',     faces, ...                 % - faces
     'vertices',  vertices, ...                 % - vertices
 'facelighting',   'phong', ...                 % - material
    'linestyle',    'none', ...                 % - edges
    'facecolor', [.7 .9 .7]);                   % - color

v = 1;                                          % Vertex-Count
x = 0; y = 0; z = -1;                           % Arrow-head Position 

for i = 0:5:720                                 % Loop over Spiral
    if i ~= 0, delete(a); end                   % - delete arrow
    
    x0 = x; y0 = y; z0 = z;                     % - save old position 
    
    x = cosd(i);                                % - define new position
    y = sind(i);
    z = 2*i/720 - 1;
    
    color = [abs(cosd(i/10)*cosd(i/10)), ...    % - update line color
             abs(cosd(i/10)*sind(i/10)), ...
             abs(sind(i/10))];
    
    vertices(3*v-2:3*v, 1:3) = [x0, y0, z0;     % - add vertices
                                 0,  0, -1; 
                                 x,  y,  z];
    v = v + 1;                                  % - increase vertex count
    set(spiral, 'vertices', vertices);          % - update spiral patch
    
    a(1) = farrow(0,0,z, x,y,z, color,2);        % - radial arrow
    a(2) = farrow(0,0,1, 0,0,z, color,2);        % - up arrow
    pause(.01);
end

pause(.8);

% (5) Flower Arrow

cla;                                            % Clear axes

set(info, ...                                   % Update text
    'string', 'Flower Arrow', ...    
    'color', [.5 .5 .5] + rand(1, 3)/2);        % - new text Color

line0 = animatedline('linewidth',        2, ... % Define animated line
                         'color',[.6 .2 .5]);

color1 = rand(1, 3);    color2 = rand(1, 3);    % Define 4 random colors
color3 = rand(1, 3);    color4 = rand(1, 3);

for i = 0:3:360                                 % Loop over Flower
    
    if i ~= 0, delete(a); end                   % - delete arrows
    
    x = cosd(2*i)*cosd(i);                      % - update position
    y = cosd(2*i)*sind(i);
    addpoints(line0, x, y, 0);                  % - add to animated line
    
    a(1) = farrow( 1,-1,-1, x,y,0, color1, 2);     % - define arrows
    a(2) = farrow(-1,-1,-1, x,y,0, color2, 2); 
    a(3) = farrow(-1, 1,-1, x,y,0, color3, 2);
    a(4) = farrow( 1, 1,-1, x,y,0, color4, 2);
    
    pause(.01);
end