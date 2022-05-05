function partitionPhantom = InitPartitionPhantom(...
    numSpins, ...
    phantomWidth, ...
    compartmentHeights, ...
    T1Vec, ...
    T2Vec, ...
    M0Vec, ...
    offsetVec)
% Creates a partition function, schematically drawn below:
%
%  phantomWidth (x)
% <--------------->
%  _______________
% |               |      /|\
% |               |       | 
% |    Comp. N    |       | 
% |               |       | 
% |               |       | 
% |_______________|       | 
% |       .       |       | 
% |       .       |       | 
% |_______________|       |
% |    Comp. 3    | _ _ _ | _ _ _ _ _ _ y = 0
% |               |       | 
% |_______________|       | phantomHeight (y)
% |               |       |
% |    Comp. 2    |       | 
% |_______________|       |
% |               |       |
% |               |       |
% |    Comp. 1    |       |
% |               |       |
% |_______________|      \|/
%
%
% Each partition has its own T1, T2, M0 and offset.
%
% Inputs:
%
% Var. Name           Units   Description
% numSpins            -       A 1x2 vector, specifying the number of spins 
%                             along phantom's width and height,
%                             respectively
% phantomWidth        mm      Width of phantom
% compartmentHeights  mm      Vector of (non-negative) heights of
%                             compartments (e.g. [100, 10, 10, 100])
% T1Vec               ms      Longitudinal relaxation in each compartment
% T2Vec               ms      Transverse relaxation in each compartment
% M0Vec               -       Equilibrium magnetization of each spin
% offsetVec           kHz     Chemical shift

% Define phantom constants
phantomHeight = sum(compartmentHeights); % mm

% Compute the center positions of each of the four compartments
numCompartments = length(compartmentHeights);
compartmentCenters = zeros(numCompartments);
compartmentCenters(1) = -phantomHeight/2 + compartmentHeights(1)/2;
for curCompartment=2:numCompartments
    compartmentCenters(curCompartment) ...
        = compartmentCenters(curCompartment-1) ...
        + compartmentHeights(curCompartment-1)/2 ...
        + compartmentHeights(curCompartment)/2;
end

% Compute number of spins in each compartment along the compartment's
% height (the y-axis)
numSpinsPerCompartment = ceil(numSpins(2)*(compartmentHeights./phantomHeight));

% Create the different compartment spins
partitionPhantom = [];
for curCompartment = 1:numCompartments
    offsetX = 0;
    offsetY = compartmentCenters(curCompartment);
    offsetZ = 0;
    % The spins are equally spaced inside a compartment, without actually
    % 'touching' the bottom and top limits. For example, a compartment of 
    % height 10 mm with 3 spins will (assuming it's centered about 0)
    % have them at -5, 0, and 5 mm (and not -10, 0, 10 mm). The reason
    % for this is to avoid "spin buildup" in regions where compartments
    % meet; i.e., avoid having two spins with the same y coordinate, one
    % from one compartment and the second from another compartment.
    compartmentHeight = compartmentHeights(curCompartment)*(1 - 1/numSpinsPerCompartment(curCompartment));
    if (numSpinsPerCompartment(curCompartment)~=0)
        spins = InitPhantomSquare(phantomWidth,...
                                  compartmentHeight,...
                                  offsetX, ...
                                  offsetY, ...
                                  offsetZ, ...
                                  numSpins(1),...
                                  numSpinsPerCompartment(curCompartment),...
                                  T1Vec(curCompartment),...
                                  T2Vec(curCompartment),...
                                  M0Vec(curCompartment),...
                                  offsetVec(curCompartment));
        partitionPhantom = [partitionPhantom, spins];
    end
end
