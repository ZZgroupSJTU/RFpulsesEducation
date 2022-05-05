% ========================================================================
% This plots the thermal polarization of spins on a bloch sphere.
% ========================================================================

N=550;
theta = rand(1,N)*pi*0.06;

r = 1;
[x,y,z] = sphere(50);
x0 = 0; y0 = 0; z0 = 0;
x = x*r + x0;
y = y*r + y0;
z = z*r + z0;

% Then you can use a surface command as Patrick suggests:



phi = rand(1,N)*2*pi;

figure 
hold on

lightGrey = 0.8*[1 1 1]; % It looks better if the lines are lighter
surface(x,y,z,'FaceColor', 'blue','EdgeColor',lightGrey, 'FaceAlpha', 0.5)
shading interp
colormap bone

hold on
for idx=1:N
    x = sin(theta(idx))*cos(phi(idx));
    y = sin(theta(idx))*sin(phi(idx));
    z = cos(theta(idx));
    line([0 x], [0 y], [0 z], 'color', [0.9 0.9 0]);
    plot3(x,y,z, 'r.', 'MarkerSize', 12);
end
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);

view(36, 12);