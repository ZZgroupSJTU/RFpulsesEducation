function spins = InitPhantomSphere(radius,...
                                   numElements,...
                                   T1,...
                                   T2,...
                                   M0,...
                                   offset)
% Description: Initializes a 3D spherical phantom of a given radius,
% with a given number of elements along each axis. 
% All the magnetization vectors start off at 
% thermal equilibrium, M=[0; 0; M0].
%
% Inputs:
%
% Var. Name   Units   Description
% radius      mm      Radius of sphere. Also serves to determine the
%                     "FOV" of the phantom: 2*radius.
% numElements -       Number of elements along the x/y/z axes.
% T1          ms      Longitudinal relaxation
% T2          ms      Transverse relaxation
% M0          a.u.    Equilibrium magnetization
% offset      kHz     Chemical shift of spins

% Create position vectors
xVec = linspace(-radius, radius, numElements);
yVec = linspace(-radius, radius, numElements);
zVec = linspace(-radius, radius, numElements);

% Create sphere
counter = 0;
for cx=1:numElements
    for cy=1:numElements
        for cz=1:numElements
            curX = xVec(cx);
            curY = yVec(cy);
            curZ = zVec(cz);
            if (sqrt(curX^2+curY^2+curZ^2)<radius)
                counter = counter + 1;
                spins(counter).r = [curX; curY; curZ];
                spins(counter).M = [0; 0; M0];
                spins(counter).cs = offset;
                spins(counter).T1 = T1;
                spins(counter).T2 = T2;
                spins(counter).M0 = M0;
                spins(counter).B1 = 1; % Scales RF
                spins(counter).B0 = 0; % Offset, in kHz
                spins(counter).RS = 1; % Receiver sensitivity
            end
        end
    end
end

