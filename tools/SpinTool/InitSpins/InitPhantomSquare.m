function spins = InitPhantomSquare(FOVX,...
                                   FOVY,...
                                   offsetX, ...
                                   offsetY, ...
                                   slicePosZ,...
                                   Nx,...
                                   Ny,...
                                   T1,...
                                   T2,...
                                   M0,...
                                   offset)
% Description: Initializes a square phantom of 
% equal magnetization along the x and y axes. 
% All the magnetization vectors start off at 
% thermal equilibrium, M=[0; 0; M0].
%
% Inputs:
%
% Var. Name   Units   Description
% FOVX        mm      Phantom size along x
% FOVY        mm      Phantom size along y
% offsetX     mm      Offset of phantom along x
% offsetY     mm      Offset of phantom along y
% slicePosZ   mm      Position of slice along z-axis
% Nx          -       Num. of spins along x
% Ny          -       Num. of spins along y
% T1          ms      Longitudinal relaxation
% T2          ms      Transverse relaxation
% M0          a.u.    Equilibrium magnetization
% offset      kHz     Chemical shift of spins


if (Nx==1)
    xx = 0;
else
    xx = linspace(-FOVX/2, FOVX/2, Nx);
end
xx = xx + offsetX;

if (Ny==1)
    yy = 0;
else
    yy = linspace(-FOVY/2, FOVY/2, Ny);
end
yy = yy + offsetY;

counter = 0;
for k=1:Nx
    for p=1:Ny
        counter = counter + 1;
        spins(counter).r = [xx(k); yy(p); slicePosZ];
        spins(counter).M = [0; 0; M0];
        spins(counter).cs = offset;
        spins(counter).T1 = T1;
        spins(counter).T2 = T2;
        spins(counter).M0 = M0;
        spins(counter).B1 = 1;
        spins(counter).B0 = 0;
        spins(counter).RS = 1;
    end;
end;