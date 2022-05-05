function partitionPhantom = InitPartitionPhantomNYU(numSpins, field, centerFreq, waterSuppressMode)
% Creates a partition phantom that matches the one found in NYU's CBI's
% Oded Gonen's laboratory. It is sketched below (w/ estimated T1s, T2s):
%
%       16 cm
% <--------------->             T1 (ms)    T2 (ms)    Offset (ppm)     M0
%  _______________
% |               |   /|\
% |               |    | 
% |     Water     |    |        1500       200        -4.7             1000
% |               |    | 
% |               |    | 
% |_______________|    | 
% |               |    | 
% |     NAA       |    |        1360       343        -2.02            2
% |_______________|    | 20 cm  
% |               |    |                      
% |    Choline    |    |        1145       248        -3.2             1
% |_______________|    | 
% |               |    |
% |   Creatine    |    |        1300       172        -3               1
% |_______________|    |
% |               |    |       
% |  myoInositol  |    |        1170       100        -3.56            0.5
% |_______________|    |       
% |               |    |
% |               |    |
% |               |    |
% |     Water     |    |        1500       200        -4.7             1000
% |               |    |
% |_______________|   \|/
%
%
% Inputs
% Variable           Units     Description
% numSpins            -        A 1x2 vector, specifying the number of spins 
%                              along phantom's width and height, respectively
% field              Tesla     The main B0 field 
% centerFreq         ppm       Where the 0 of the spectrum is
% waterSuppressMode  -         A string:
%                              'none'    - No water suppression (default)
%                              'partial' - Some "water suppression" - sets M0 of water to 5
%                              'full'    - Full "water suppression" - sets M0 of water to 0

gyromagneticRatio = GetGyromagneticRatio('1h'); % MHz/T = kHz/mT

% Phantom dimensions. Each metabolite's compartment is 1 cm thick
phantomWidth = 160; % mm
compartmentHeights = [80 10 10 10 10 80];

%             Water   NAA   Cho  Cre   myI     Water
offsetVec = [ -4.7   -2.02 -3.2  -3   -3.56   -4.7];  % in ppm
offsetVec = offsetVec - centerFreq;
offsetVec = offsetVec*(gyromagneticRatio*field*0.001); % in kHz


% Compartment properties
%          Water   NAA  Cho   Cre   myI   Water
T1Vec = [  1500,  1360, 1145, 1300, 1170, 1500];
T2Vec = [   200,   343,  248,  172,  100,  200];
switch(waterSuppressMode)
    case 'none'
        M0Vec = [1000 2 1 1 0.5 1000];
    case 'partial'
        M0Vec = [   5 2 1 1 0.5    5];
    case 'full'
        M0Vec = [   0 2 1 1 0.5    0];
    otherwise  % Default: none
        M0Vec = [1000 2 1 1 0.5 1000];
end
        
partitionPhantom = InitPartitionPhantom(...
    numSpins, ...
    phantomWidth, ...
    compartmentHeights, ...
    T1Vec, ...
    T2Vec, ...
    M0Vec, ...
    offsetVec);