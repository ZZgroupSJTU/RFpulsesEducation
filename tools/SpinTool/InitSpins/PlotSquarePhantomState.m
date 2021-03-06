function PlotSquarePhantomState(FOVX, FOVY, Nx, Ny, spins)

% Description: Plots Mz as a function of x and
% y for a given square phantom. The assumption
% is that the phantom is arranged "columns 
% first";, e.g., 
%
%  ______________________
% | 1  |  4  |  7  | 10  |
% |____|_____|_____|_____|
% | 2  |  5  |  8  | 11  |
% |____|_____|_____|_____|
% | 3  |  6  |  9  | 12  |
% |____|_____|_____|_____|
%
% Inputs:
%
% Var. Name   Units   Description
% FOVX        mm      Phantom size along x
% FOVY        mm      Phantom size along y
% Nx          -       Num. of spins along x
% Ny          -       Num. of spins along y

xx = linspace(-FOVX/2, FOVX/2, Nx);
yy = linspace(-FOVY/2, FOVY/2, Ny);

Mz = zeros(Nx, Ny);
counter = 0;
for k=1:Nx
    for p=1:Ny
        counter = counter + 1;
        Mz(k,p) = spins(counter).M(3);
    end;
end;

figure
imagesc(xx, yy, Mz, [-1 1]);
colorbar
title('Mz');
xlabel('x (mm)');
ylabel('y (mm)');