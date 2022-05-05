function area = voronoidens(kx,ky);
%
% function area = voronoidens(kx,ky);
%
% input:  kx, ky are k-space trajectories
% output: area of cells for each point 
%           (if point doesn't have neighbors the area is NaN)

[row,column] = size(kx);

% uncomment these to plot voronoi diagram
%[vx, vy] = voronoi(kx,ky);
%plot(kx,ky,'r.',vx,vy,'b-'); axis equal

kxy = [kx(:),ky(:)];
% returns vertices and cells of voronoi diagram
[V,C] = voronoin(kxy); 
area = [];
for j = 1:length(kxy)
  x = V(C{j},1); y = V(C{j},2); lxy = length(x);
  A = abs(sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:))));
  area = [area A];
end

area = reshape(area,row,column);

