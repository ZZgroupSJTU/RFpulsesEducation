function imgMag = SimulateVoxelProfile(varargin)
% Simulates the voxel profile of a localization sequence
%   [voxelProfileRE, imgProfileIM] = SimulateVoxelProfile('name', value, ... ) 
%   Returns two VDIImages, with the real and imaginary parts of the voxel 
%   profile, as deduced using Bloch simulations. Parameter name-value pairs
%   allow tailoring of function to specific sequences, pulses, etc:
%     pulses      A cell array of the selective pulses. The number of
%                 elements is determined by the sequence. For example, a
%                 PRESS will require 3 pulses. The pulses should already
%                 have their gradients calibrated. THESE MUST BE SUPPLIED.
%     seq         Default: 'PRESS'. Current options: 'PRESS', 'STEAM', 
%                 'STRESS'.
%     TE          Vector of inter-pulse delays, with as many elements as
%                 there are pulses. Default: sequence dependent. For PRESS,
%                 STEAM and STRESS, the minimum TE attainable is used,
%                 as calculated by the function CalcMinEchoTime.
%     FOV         A 1x3 vector, in mm, of field of view of the resulting
%                 images. Default: [32 32 32].
%     res         A 1x3 vector of number of voxels in resulting images.
%                 Default: [32 32 32].
%     offset      Default: [0 0 0]. Denotes the center of the resulting
%                 images (in mm). 
%     orientation A 3x3 matrix defining the orientation of the produced
%                 images. Default: eye(3) (3x3 identity matrix).
%     img         A VDIImage which is used as a template. If provided,
%                 the size, offset and resolution are taken directly from
%                 the image, and "FOV", "res" and "offset" above are 
%                 ignored.
%     B1          B1 scaling (to model inhomogeneity). Default: 1. 
%     T1          T1 of spins. Default: inf. Could be a number (in ms),
%     T2          T2 of spins. Default: inf. Could be a number (in ms),
%
%   NOTE: If B1, T1 or T2 are provided as VDIImages, they should match the
%   
p = inputParser;
p.addParameter('pulses', [], @(x) iscell(x));
p.addParameter('seq', 'PRESS', @(x) ismember(lower(x), {'press', 'steam', 'stress'}));
p.addParameter('FOV', [32; 32; 32], @(x) isvector(x));
p.addParameter('orientation', eye(3), @(x) ismatrix(x));
p.addParameter('res', [32; 32; 32], @(x) isvector(x));
p.addParameter('offset', [0;0;0], @(x) isvector(x));
p.addParameter('img', [], @(x) isa(x, 'VDIImage'));
p.addParameter('B1', 1, @(x) isnumeric(x));
p.addParameter('T1', inf, @(x) isnumeric(x));
p.addParameter('T2', inf, @(x) isnumeric(x));
p.addParameter('TE', [], @(x) isempty(x) || isvector(x));
p.CaseSensitive = false;
parse(p, varargin{:});
inputParams = p.Results;

if ~isempty(inputParams.img)
    offset = img.position;
    orientation = img.orientation;
    res = img.matrixRes;
    FOV = img.size;
else
    offset = inputParams.offset;
    orientation = inputParams.orientation;
    res = inputParams.res;
    FOV = inputParams.FOV;
    if isrow(offset), offset = offset.'; end
    if isrow(res), res = res.'; end
    if isrow(FOV), FOV = FOV.'; end
end

pulseDurations = [inputParams.pulses{1}.tp, inputParams.pulses{2}.tp, inputParams.pulses{3}.tp];
switch lower(inputParams.seq)
    case 'press'
        numPulses = 3;
        [~, minTEPulseSpacings] = CalcMinEchoTime('sequence', 'press', 'pulseDurations', pulseDurations);
    case 'steam'
        numPulses = 3;
        % [~, minTEPulseSpacings] = CalcMinEchoTime('sequence', 'steam', 'pulseDurations', pulseDurations);
        error('Sequence not yet implemented.');
    case 'stress'
        numPulses = 3;
        % [~, minTEPulseSpacings] = CalcMinEchoTime('sequence', 'stress', 'pulseDurations', pulseDurations);
        error('Sequence not yet implemented.');
    otherwise
        numPulses = 3;
        % [~, minTEPulseSpacings] = CalcMinEchoTime('sequence', 'press', 'pulseDurations', pulseDurations);
        error('Unrecognized sequence string.');
end

if isempty(inputParams.TE)
    TE = minTEPulseSpacings;
else
    TE = inputParams.TE;
end

if numel(TE)~=numPulses
    error('Number of echo times (%d) does not match number of expected pulses (%d) for sequence (%s)', numel(TE), numPulses, inputParams.seq);
end

interPulseDelay = zeros(1, numPulses);
for idxPulse=1:numPulses-2
    interPulseDelay(idxPulse) = TE(idxPulse) - pulseDurations(idxPulse)/2 - pulseDurations(idxPulse+1)/2;
end
interPulseDelay(numPulses) = pulseDurations(numPulses)/3;

% =========================================================================
% Create spins
% =========================================================================

imgMz = VDIImage(FOV, orientation, offset, ones(res.'));

% Populate spin structure
spins = InitSpins3D('numSpins', res.', 'sampleSize', FOV.', 'sampleCenter', offset.', 'T1', inputParams.T1, 'T2', inputParams.T2, 'B1',inputParams.B1);

% =========================================================================
% Simulate voxel profile
% =========================================================================

switch lower(inputParams.seq)
    case 'press'
        refocMomentX = -inputParams.pulses{1}.tp/2*inputParams.pulses{1}.Gx(1)*1000;
        refocMomentY = -inputParams.pulses{1}.tp/2*inputParams.pulses{1}.Gy(1)*1000;
        refocMomentZ = -inputParams.pulses{1}.tp/2*inputParams.pulses{1}.Gz(1)*1000;
        seq = {{'pulse', inputParams.pulses{1}}, {'purgemoment',  refocMomentX, refocMomentY, refocMomentZ, interPulseDelay(1)+1e-6}, ...
               {'pulse', inputParams.pulses{2}, [0 90 180 270], [1 1 1 1]/4, 'inverted'}, {'delay',  interPulseDelay(2)}, ...
               {'pulse', inputParams.pulses{3}, [0 90 180 270], [1 1 1 1]/4, 'inverted'}, {'delay',  interPulseDelay(3)}};
        spinsOut = ApplySequence3D(spins, seq);   
    case 'steam'
        error('Sequence not yet implemented.');
    case 'stress'
        error('Sequence not yet implemented.');
    otherwise
        error('Unrecognized sequence string.');
end

imgMag{1} = imgMz.Copy;  imgMag{1}.data = spinsOut.M(:,:,:,1);  imgMag{1}.name = 'Mx';    imgMag{1}.clim = [-1 1];
imgMag{2} = imgMz.Copy;  imgMag{2}.data = spinsOut.M(:,:,:,2);  imgMag{2}.name = 'My';    imgMag{2}.clim = [-1 1];
imgMag{3} = imgMz.Copy;  imgMag{3}.data = spinsOut.M(:,:,:,3);  imgMag{3}.name = 'Mz';    imgMag{3}.clim = [-1 1];
