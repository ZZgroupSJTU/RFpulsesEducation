function spins = InitSpins3D(varargin)
% Initializes a 3D rectilinear spin structure for simulations
%   spins = InitSpins3D('param', value, ... ) Possible parameter names 
%   are:
%     numSpins      1x3 array of number of spins along x, y, z. 
%                   Default: 11*11*11
%     sampleSize    1x3 array of size of spin array along x, y, z axes,
%                   in mm. Default: 10*10*10. 
%     sampleCenter  1x3 array of center of 3D array in (x,y,z) space, in mm.
%                   Default: [0 0 0].
%     T1            T1 value, shared by all spins, in ms. Default: 1 sec.
%     T2            T2 value, shared by all spins, in ms. Default: 0.1 sec.
%     cs            Chemical shift offset, shared by all spins, in kHz. 
%                   Default: 0.
%     B1            B1 scaling, shared by all spins. Default: 1.
%     M0            Equilibrium z-magnetization, shared by all spins. 
%                   Default: 1.
%     initMag       Initial magnetization. This can either by a single number,
%                   which will be assumed to apply to the Mz for all spins;
%                   a 3x1 or 1x3 vector, which will be assumed to be the
%                   magnetization vector which applies to all spins; or a
%                   full Nx*Ny*Nz*3 matrix, where Nx*Ny*Nz matches 
%                   the numSpins input. Default: initMag = 1.

p = inputParser;
p.addParameter('numSpins', [10 10 10], @(x) isnumeric(x) && numel(x)==3); 
p.addParameter('sampleSize', [10 10 10], @(x) isnumeric(x) && numel(x)==3); 
p.addParameter('sampleCenter', [0 0 0], @(x) isnumeric(x) && numel(x)==3); 
p.addParameter('T1', 1000, @(x) isnumeric(x) && numel(x)==1); 
p.addParameter('T2', 100, @(x) isnumeric(x) && numel(x)==1); 
p.addParameter('cs', 0, @(x) isnumeric(x) && numel(x)==1); 
p.addParameter('B1', 1, @(x) isnumeric(x) && numel(x)==1); 
p.addParameter('M0', 1, @(x) isnumeric(x) && numel(x)==1); 
p.addParameter('initMag', 1, @(x) isnumeric(x)); 
p.parse(varargin{:});

spins.T1 = p.Results.T1;
spins.T2 = p.Results.T2;
spins.cs = p.Results.cs;
spins.B1 = p.Results.B1;
spins.M0 = p.Results.M0;

sampleSize = p.Results.sampleSize;
sampleCenter = p.Results.sampleCenter;
numSpins = p.Results.numSpins;
spins.xVec = linspace(sampleCenter(1) - sampleSize(1)/2, sampleCenter(1) + sampleSize(1)/2, numSpins(1));
spins.yVec = linspace(sampleCenter(2) - sampleSize(2)/2, sampleCenter(2) + sampleSize(2)/2, numSpins(2));
spins.zVec = linspace(sampleCenter(3) - sampleSize(3)/2, sampleCenter(3) + sampleSize(3)/2, numSpins(3));

% Let's check how the initial magnetization is supplied and act accordingly.
s = size(p.Results.initMag);
initMag = p.Results.initMag;
switch numel(p.Results.initMag)
    case 1
        spins.M = zeros([numSpins 3]);
        spins.M(:,:,:,3) = initMag;
    case 3
        spins.M = zeros([numSpins 3]);
        spins.M(:,:,:,1) = initMag(1);
        spins.M(:,:,:,2) = initMag(2);
        spins.M(:,:,:,3) = initMag(3);
    otherwise
        if (~isequal(s(1:3), numSpins(1:3))) || (s(4)~=3)
            error('The initial magnetization, if supplied as an array, should match the number of spins used.');
        end
        spins.M = initMag;
end


