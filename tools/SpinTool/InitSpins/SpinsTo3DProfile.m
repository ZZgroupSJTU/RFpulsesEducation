function [xAxis, yAxis, zAxis, magnetization] = SpinsTo3DProfile(spins, magIndex)
% SpinsTo3DProfile  Convert a spin structure to a 3D profile 
%
%   [xAxis, yAxis, zAxis, magnetization] = SpinsTo3DProfile(spins)
% 
% Inputs 
%
% Name          Type         Units   Description      
% spins         -            -       Input spin structure, 
%                                    with relaxation
% magIndex      Integer      -       = 1, 2, 3 for Mx, My or Mz
%
% Outputs
%
% Name          Type         Units   Description      
% x/y/zAxis     Vectors      mm      Axes for plotting data
% magnetization 3D array     a.u.    3D data containing Mz

numSpins = numel(spins);
posMat = zeros(3,numSpins);
data = zeros(1,numSpins);
for idx=1:numSpins
    posMat(:,idx) = spins(idx).r;
    data(idx) = spins(idx).M(magIndex);
end

% These axes are monotonically increasing, but not necessarily equispaced
xPos = posMat(1,:)';
yPos = posMat(2,:)';
zPos = posMat(3,:)';
F = TriScatteredInterp(xPos, yPos, zPos, data');

xAxis = unique(posMat(1,:))';
yAxis = unique(posMat(2,:))';
zAxis = unique(posMat(3,:))';
[xMat, yMat, zMat] = meshgrid(xAxis, yAxis, zAxis);

magnetization = F(xMat, yMat, zMat);
