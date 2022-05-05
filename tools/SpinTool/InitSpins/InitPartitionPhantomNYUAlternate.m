function partitionPhantom = InitPartitionPhantomNYUAlternate(numSpins, field)
% Creates a partition phantom that matches the one found in NYU's CBI's
% Oded Gonen's laboratory. It is sketched below. It differs from the actual
% partition phantom by the following attributes:
% 1. T1 and T2 of all metabolites are set to 10000 milliseconds, to
%    eliminate relaxation effects.
% 2. Water amplitude is set to 0.
% 3. All peaks are spaced 1 ppm apart, at -1.5, -0.5, 0.5, 1.5 ppm.
% 4. All peaks have the same amplitude.
%
%       16 cm
% <--------------->            
%  _______________
% |               |   /|\
% |               |    | 
% |     Water     |    |       
% |               |    | 
% |               |    | 
% |_______________|    | 
% |               |    | 
% |     NAA       |    |       
% |_______________|    | 20 cm  
% |               |    |                      
% |    Choline    |    |       
% |_______________|    | 
% |               |    |
% |   Creatine    |    |       
% |_______________|    |
% |               |    |       
% |  myoInositol  |    |       
% |_______________|    |       
% |               |    |
% |               |    |
% |               |    |
% |     Water     |    |       
% |               |    |
% |_______________|   \|/
%
%
% Inputs
% Variable           Units     Description
% numSpins            -        A 1x2 vector, specifying the number of spins 
%                              along phantom's width and height, respectively
% field              Tesla     The main B0 field 

gyromagneticRatio = 42.57; % MHz/T = kHz/mT

% Phantom dimensions. Each metabolite's compartment is 1 cm thick
phantomWidth = 160; % mm
compartmentHeights = [80 10 10 10 10 80];

%             Water   NAA   Cho  Cre   myI     Water
offsetVec = [   0    -1.5  -0.5  0.5   1.5       0];  % in ppm
offsetVec = offsetVec*(gyromagneticRatio*field*0.001); % in kHz


% Compartment properties
%          Water   NAA  Cho   Cre   myI   Water
T1Vec = [  100000,  100000, 100000, 100000, 100000, 100000];
T2Vec = [  100000,  100000, 100000, 100000, 100000, 100000];
M0Vec = [   0   ,      1  ,   1  ,     1  ,   1  ,     0  ];

partitionPhantom = InitPartitionPhantom(...
    numSpins, ...
    phantomWidth, ...
    compartmentHeights, ...
    T1Vec, ...
    T2Vec, ...
    M0Vec, ...
    offsetVec);