function spins = InitPhantomCircle(radius,...
                                     posAlongNormal, ...
                                     numElements,...
                                     T1,...
                                     T2,...
                                     M0,...
                                     offset, ...
                                     normal)
% Description: Initializes a 2D circular phantom (in the XY plane) 
% of a given radius, with a given number of elements along each axis. 
% All the magnetization vectors start off at thermal eq.: M=[0; 0; M0].
%
% Inputs:
%
% Var. Name         Units          Description
% radius            mm             Radius of sphere. Also serves to determine the
%                                  "FOV" of the phantom: 2*radius.
% posAlongNormal    mm             z ordinate of spins (same for all spins)
% numElements       -              Number of elements along the x/y/z axes.
% T1                ms             Longitudinal relaxation
% T2                ms             Transverse relaxation
% M0                a.u.           Equilibrium magnetization
% offset            kHz            Chemical shift of spins
% normal            'x', 'y', 'z'  Normal to circle's plane

% Create position vectors
xVec = linspace(-radius, radius, numElements);
yVec = linspace(-radius, radius, numElements);
zVec = linspace(-radius, radius, numElements);

% Create sphere
counter = 0;
for idxA=1:numElements
    for idxB=1:numElements
        switch (normal)
            case 'x'
                curX = posAlongNormal;
                curY = yVec(idxA);
                curZ = zVec(idxB);
                isInside = ((curY^2 + curZ^2) < radius^2);
            case 'y'
                curX = xVec(idxA);
                curY = posAlongNormal;
                curZ = zVec(idxB);
                isInside = ((curZ^2 + curX^2) < radius^2);
            case 'z'
                curX = xVec(idxA);
                curY = yVec(idxB);
                curZ = posAlongNormal;
                isInside = ((curX^2 + curY^2) < radius^2);
        end
        if isInside
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

