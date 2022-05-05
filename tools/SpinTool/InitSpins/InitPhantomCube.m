function spins = InitPhantomCube(FOV, offset, numSpins, T1, T2, M0, chemShift, B1)
% Description: Initializes a square phantom of 
% equal magnetization along the x and y axes. 
% All the magnetization vectors start off at 
% thermal equilibrium, M=[0; 0; M0].
%
% Inputs:
%
% Var. Name           Units   Description
% FOVX, FOVY, FOVZ    mm      Phantom size along x, y or z
% offsetX, offsetY,           
%          offsetZ    mm      Offset of phantom along x, y or z
% Nx, Ny, Nz          -       Num. of spins along x, y or z
% T1                  ms      Longitudinal relaxation
% T2                  ms      Transverse relaxation
% M0                  a.u.    Equilibrium magnetization
% offset              kHz     Chemical shift of spins
% B1                  -       B1 scaling factor

FOVX = FOV(1);
FOVY = FOV(2);
FOVZ = FOV(3);
 
Nx = numSpins(1);
Ny = numSpins(2);
Nz = numSpins(3);

offsetX = offset(1);
offsetY = offset(2);
offsetZ = offset(3);

dx = FOVX/Nx;
dy = FOVY/Ny;
dz = FOVZ/Nz;

if (Nx==1)
    xx = 0;
else
    xx = linspace(-FOVX/2+dx/2, FOVX/2-dx/2, Nx);
end
xx = xx + offsetX;

if (Ny==1)
    yy = 0;
else
    yy = linspace(-FOVY/2+dy/2, FOVY/2-dy/2, Ny);
end
yy = yy + offsetY;

if (Nz==1)
    zz = 0;
else
    zz = linspace(-FOVZ/2+dz/2, FOVZ/2-dz/2, Nz);
end
zz = zz + offsetZ;


counter = 0;
numSpins = Nx*Ny*Nz;
spins(numSpins).r=[0;0;0];
for k=1:Nx
    for p=1:Ny
        for m=1:Nz
            counter = counter + 1;
            spins(counter).r = [xx(k); yy(p); zz(m)];
            spins(counter).M = [0; 0; M0];
            spins(counter).cs = chemShift;
            spins(counter).T1 = T1;
            spins(counter).T2 = T2;
            spins(counter).M0 = M0;
            spins(counter).B1 = B1; % Scales RF
            spins(counter).B0 = 0; % Offset, in kHz
            spins(counter).RS  = 1; % Receiver sensitivity
        end
    end
end