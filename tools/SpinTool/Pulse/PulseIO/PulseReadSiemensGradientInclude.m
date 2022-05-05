function gradShape = PulseReadSiemensGradientInclude(filename)
% SYNTAX: 
%
%     gradShape = PulseReadSiemensGradientInclude(filename)
%
% Reads a (normalized) gradient shape from a .h file.
%
% Inputs:
%
% Variable Name   Units    Description
% filename        -        Full filename, including directory. 
%
% Outputs:
%
% Variable Name   Units    Description
% gradShape       -        The normalized gradient shape.

filename = RemoveExtensionFromFilename(filename);
filename = [filename, '.h'];


% -----------------------------------------------------------------------
% Load pulse data
% ----------------------------------------------------------------------

% Load data
fid=fopen(filename, 'r'); 
pulseDataStr = fread(fid, inf, 'uint8=>char')';
fclose(fid);

odds = strfind(pulseDataStr, '(');
evens = strfind(pulseDataStr, ')');
numPoints = numel(odds);



gradShape = zeros(1, numPoints);

for idx=1:numPoints
    gradShape(idx) = str2double(pulseDataStr(odds(idx)+1:evens(idx)-1));
end

